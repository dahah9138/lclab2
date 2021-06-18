#include <lclab2.h>

#include "Widget.h"

#define USE_PLANE 1


using namespace Magnum;
using namespace Math::Literals;

/*
    TODO: Make sheet color draw dynamically

    Next: Maybe consider making a POM imaging lib for lclab2 (make with CPU first)
    - Needs Jones Matrix (I have matrices through Eigen)
    - Needs Lamp specifications
    - LC Optical properties
*/

class Sandbox : public LC::Application
{
public:

	explicit Sandbox(const Arguments& arguments);

	~Sandbox();

private:
	virtual void drawEvent() override;

    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;
    void textInputEvent(TextInputEvent& event) override;
    void keyPressEvent(KeyEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;
    void updateColor();

    // Tested geometries
    LC::SphereArray _grid;
    LC::Torus _sheet;
    LC::DynamicColorSheet _dsheet;
    // Torus with PhongGL shader
    LC::NormalTorus _sheetNormal;

    /*
        GUI Widget
        *** Make sure to keep decoupled from actual simulation core functionality
    */
    Widget _widget;

    // Solver
    LC::Solver* _solver;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("Sandbox Application")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable) 
                                                              } {
    /* Setup the GUI */
    setupGUI();

    /* Setup window and parameters */
    enableDepthTest();
    enableFaceCulling();

    /* Loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
    

    /* Setup camera */
    setupCamera(0.9f);

    
    //_sheet.Init();
    //_sheetNormal.Init();


    _solver = new LC::FrankOseen::ElasticOnly::FOFDSolver;

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;
    using T4 = LC::FrankOseen::ElasticOnly::FOFDSolver::Tensor4;

    /* Setup data */
    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    data->voxels[0] = 32;
    data->voxels[1] = 32;
    data->voxels[2] = 32;

    // hard
    data->bc[0] = 0;
    data->bc[1] = 0;
    data->bc[2] = 0;

    // toron stable dimensions
    data->cell_dims[0] = 0.35;
    data->cell_dims[1] = 0.35;
    data->cell_dims[2] = 0.35;

    data->k11 = LC::FrankOseen::ElasticConstants::_5CB("k11");
    data->k22 = LC::FrankOseen::ElasticConstants::_5CB("k22");
    data->k33 = LC::FrankOseen::ElasticConstants::_5CB("k33");

    // Generate a toron
    // Make a director modifying config that takes a function pointer to a lambda that modifies that director

    auto toron = [](T4 &n, int i, int j, int k, int *voxels) {
    
        int d[3] = { voxels[0] / 4, voxels[1] / 4, voxels[2] / 4 };

        if (abs(k - voxels[2] / 2) < d[2] && abs(i - voxels[0] / 2) < d[0] && abs(j - voxels[1] / 2) < d[1]) {
            n(i, j, k, 2) = -1.0;
        }
        else {
            n(i, j, k, 2) = 1.0;
        }

        n(i, j, k, 0) = 0.0;
        n(i, j, k, 1) = 0.0;
        
    };

    data->config = toron;

    _solver->Init();

    /* Setup spheres */
    _grid.NX = data->voxels[0];
    _grid.NY = data->voxels[1];
    _grid.CX = data->cell_dims[0];
    _grid.CY = data->cell_dims[1];

    _grid.Init();

    _dsheet.NX = data->voxels[0];
    _dsheet.NY = data->voxels[1];
    _dsheet.CX = data->cell_dims[0];
    _dsheet.CY = data->cell_dims[1];
    _dsheet.Init();


    // Colors
    updateColor();



    LC_INFO("Created Sandbox!");
}

Sandbox::~Sandbox() {
	LC_INFO("Destroying Sandbox.");

    delete _solver;
}

/*
    Main simulation loop
*/
void Sandbox::drawEvent()
{

    GL::defaultFramebuffer.clear(
        GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    _imgui.newFrame();

    /* 1. Show a simple window.
       Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets appear in
       a window called "Debug" automatically */
    {
        ImGui::Text("Hello, world!");
        ImGui::SliderFloat("Float", &_widget.floatValue, 0.0f, 1.0f);
        if (ImGui::ColorEdit3("Clear Color", _widget.clearColor.data()))
            GL::Renderer::setClearColor(_widget.clearColor);
        if (ImGui::Button("Test Window"))
            _widget.showDemoWindow ^= true;
        if (ImGui::Button("Another Window"))
            _widget.showAnotherWindow ^= true;

        // Set Cycle
        ImGui::InputInt("Cycle", &_widget.cycle);

        // Relaxation rate
        {
            using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;

            Dataset* data = (Dataset*)(_solver->GetDataPtr());
            Float relaxRate = data->rate;
            ImGui::InputFloat("Relax rate", &relaxRate);
            data->rate = relaxRate;
        }

        // Pressed the relax button
        _widget.relax = ImGui::Button("Relax");
        ImGui::SameLine();
        _widget.print = ImGui::Button("Print");
        
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
            1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
    }
    /* 2. Show another simple window, now using an explicit Begin/End pair */
    if (_widget.showAnotherWindow) {
        ImGui::SetNextWindowSize(ImVec2(500, 100), ImGuiCond_FirstUseEver);
        ImGui::Begin("Another Window", &_widget.showAnotherWindow);
        ImGui::Text("Hello");
        ImGui::End();
    }

    /* 3. Show the ImGui demo window. Most of the sample code is in
       ImGui::ShowDemoWindow() */
    if (_widget.showDemoWindow) {
        ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
        ImGui::ShowDemoWindow();
    }

    /* Update application cursor */
    _imgui.updateApplicationCursor(*this);


    /* Update camera */
    const bool moving = _arcballCamera->updateTransformation();

    /* Reset state. Only needed if you want to draw something else with
        different state after. */

    polyRenderer();

#if USE_PLANE
    _dsheet.Draw(_arcballCamera, _projectionMatrix);
#else
    _grid.Draw(_arcballCamera, _projectionMatrix);
#endif
    //_sheet.Draw(_arcballCamera, _projectionMatrix);
    //_sheetNormal.Draw(_arcballCamera, _projectionMatrix);

    // Make sure to draw gui last, otherwise the graphics will write over the GUI
    {
        /* Set appropriate states. If you only draw ImGui, it is sufficient to
           just enable blending and scissor test in the constructor. */
        guiRenderer();

        _imgui.drawFrame();
    }

    swapBuffers();

    if (_widget.relax) {

        _solver->Relax(_widget.cycle);

        // Update sphere colors here
        // For now just xy cross section
        updateColor();
    }
    
    
    if (_widget.print) {
        _solver->Print();
    }


    if (moving || _ioUpdate) redraw();
}

// Todo: inherit all these events through Application
void Sandbox::mousePressEvent(MouseEvent& event) {

    //for (std::size_t i = 0; i < _grid.spherePositions.size(); ++i) {
    //    const Vector3 tmpCol = Vector3(std::rand(), std::rand(), std::rand()) /
    //        Float(RAND_MAX);
    //    _grid.sphereInstanceData[i].color = tmpCol;
    //}

    _imgui.handleMousePressEvent(event);

    if (!_io->WantCaptureMouse) {
        _arcballCamera->initTransformation(event.position());
    }
    else _ioUpdate = true;

}

void Sandbox::textInputEvent(TextInputEvent& event) {
    if (_imgui.handleTextInputEvent(event)) _ioUpdate = true;
}

void Sandbox::mouseReleaseEvent(MouseEvent& event) {

    if (_imgui.handleMouseReleaseEvent(event)) _ioUpdate = true;
}

void Sandbox::keyPressEvent(KeyEvent& event) {
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::keyReleaseEvent(KeyEvent& event) {
    if (_imgui.handleKeyReleaseEvent(event)) _ioUpdate = true;
}

void Sandbox::updateColor() {

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;

    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    // Colors
    std::size_t slice = data->voxels[0];
    std::size_t cross_slice = data->voxels[1] * slice;
    std::size_t volslice = data->voxels[2] * cross_slice;

    // Convenient indexing lambdas
    auto global_matlab_idx = [volslice, cross_slice, slice](int i, int j, int k, int l) {

        return volslice * l + cross_slice * k + slice * j + i;
    };

    auto cross_idx = [slice](int i, int j) {
        return slice * j + i;
    };

    LC::scalar theta, phi, nx, ny, nz;
    for (int i = 0; i < data->voxels[0]; i++) {
        for (int j = 0; j < data->voxels[1]; j++) {

            nx = data->directors[global_matlab_idx(i, j, data->voxels[2] / 2, 0)];
            ny = data->directors[global_matlab_idx(i, j, data->voxels[2] / 2, 1)];
            nz = data->directors[global_matlab_idx(i, j, data->voxels[2] / 2, 2)];
            // Compute theta and phi

            theta = acos(nz);
            phi = M_PI + atan2(ny, nx);

            // Compute color

            Color3 hsv;
            hsv[0] = phi / (2.0f * M_PI);
            hsv[1] = theta / M_PI * 2.0f;
            hsv[2] = 2.0f - theta / M_PI * 2.0f;

            if (hsv[0] > 1.0f) hsv[0] = 1.0f;
            if (hsv[1] > 1.0f) hsv[1] = 1.0f;
            if (hsv[2] > 1.0f) hsv[2] = 1.0f;

            _grid.sphereInstanceData[cross_idx(i, j)].color = Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] });
            _dsheet.data[cross_idx(i, j)].color = Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] });
        }
    }

    // Update sheet
    _dsheet.sheetBuffer.setData(_dsheet.data, GL::BufferUsage::DynamicDraw);
}

LC::Application* LC::createApplication(int argc, char **argv) {

	return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}