#include <lclab2.h>

#include "Widget.h"

using namespace Magnum;
using namespace Math::Literals;

/*
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

    // Tested geometries
    LC::SphereArray _grid;
    LC::Torus _sheet;
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

    /* Setup spheres */
    //_grid.Init();
    //_sheet.Init();
    _sheetNormal.Init();


    _solver = new LC::FrankOseen::ElasticOnly::FOFDSolver;

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;

    /* Setup data */
    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    data->voxels[0] = 10;
    data->voxels[1] = 10;
    data->voxels[2] = 10;

    data->cell_dims[0] = 1.0;
    data->cell_dims[1] = 1.0;
    data->cell_dims[2] = 1.0;

    data->k11 = LC::FrankOseen::ElasticConstants::_5CB("k11");
    data->k22 = LC::FrankOseen::ElasticConstants::_5CB("k22");
    data->k33 = LC::FrankOseen::ElasticConstants::_5CB("k33");

    _solver->Init();


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

        // Pressed the relax button
        _widget.relax = ImGui::Button("Relax");
        
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

    //_grid.Draw(_arcballCamera, _projectionMatrix);
    //_sheet.Draw(_arcballCamera, _projectionMatrix);
    _sheetNormal.Draw(_arcballCamera, _projectionMatrix);

    {
        /* Set appropriate states. If you only draw ImGui, it is sufficient to
           just enable blending and scissor test in the constructor. */
        guiRenderer();

        _imgui.drawFrame();
    }


    swapBuffers();


    if (_widget.relax) {
        _solver->Relax(1);
    }


    // TODO: Add a button with imgui to start the relax
    bool pressedRelaxAndRelaxFinished = false;

    if (moving || _ioUpdate) redraw();
    //redraw();
}

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
    if (_imgui.handleTextInputEvent(event)) return;
}

void Sandbox::mouseReleaseEvent(MouseEvent& event) {

    if (_imgui.handleMouseReleaseEvent(event)) return;
}

void Sandbox::keyPressEvent(KeyEvent& event) {
    if (_imgui.handleKeyPressEvent(event)) return;
}

void Sandbox::keyReleaseEvent(KeyEvent& event) {
    if (_imgui.handleKeyReleaseEvent(event)) return;
}



LC::Application* LC::createApplication(int argc, char **argv) {

	return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}