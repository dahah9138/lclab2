#include "DropDownMenu.h"
#include "Widget.h"
#include "Boxes.h"
#include "Simulation.h"
#include "ZProfile.h"

using namespace Magnum;
using namespace Math::Literals;
using AppSolver = LC::FrankOseen::Electric::FOFDSolver;
using Dataset = AppSolver::Dataset;
using Geometry = LC::EllipsoidArray;


class Sandbox : public LC::Application
{
public:
    explicit Sandbox(const Arguments& arguments);
    ~Sandbox();

    void keyPressEvent(KeyEvent& event) override;

    void updateBoxes(bool first = false);
    void updateFilm();
    void updateFilmGraphics();

    void drawBoxes();
    void drawFilm();
    void drawProfile();
    void mainWindow();
    void updateCamera();
    void plotFilmData();
    Dataset* dataset();

private:
    virtual void drawEvent() override;

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader = Shaders::VertexColorGL3D{};
    SceneGraph::DrawableGroup3D _transparentDrawables;

    Boxes _boxes;
    Widget _widget;

    std::unique_ptr<Simulation> _sim;
    Containers::Optional<Geometry> _film;
    Containers::Optional<Zprofile> _zprofile;
    std::vector<unsigned int> _surfaceNodes;
    unsigned int _numSurfaceNodes;

    LC::Imaging::UniformGrid::POM _pomImager;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("Film Response Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable) } {
    Utility::Arguments args;
    args.addOption('q', "topological-charge", "1")
        .setHelp("topological-charge", "topological charge", "Q")
        .addSkippedPrefix("magnum")
        .parse(arguments.argc, arguments.argv);


    /* Setup the GUI */
    setupGUI();

    /* Setup window and parameters */
    enableDepthTest();
    enableFaceCulling();

    /* Loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);

    setupCamera(0.1f, CameraType::ArcBall);

    // Will be useful when initializing data from a file
    _solver = std::make_unique<AppSolver>();

    // Load a director file
    _header.readFile = "D:\\dev\\lclab2\\data\\simulations\\homeotropic\\smallest_toron.lmt";
    _solver->Import(_header);
    LC_INFO("Imported file {0}", _header.readFile.c_str());

    _sim = std::make_unique<Simulation>();

    // Populate sim with data
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    // Cell dims stored as d/p
    auto celldims = data->cell_dims;
    LC::scalar cvol = 1.;
    LC::scalar pitch = _widget.pitch * 1e-6; // m
    for (int d = 0; d < 3; d++) {
        celldims[d] *= pitch;
        cvol *= celldims[d] / (data->voxels[d] - 1);
    }


    // Update element masses
    // Note a 10 um x 10 um x 10 um element
    // has a mass of 10^-12 kg

    //_sim->Parameters().bulkMass *= cvol / 10e-15;
    //_sim->Parameters().surfaceMass *= cvol / 10e-15;
    //_sim->Parameters().edgeMass *= cvol / 10e-15;
    //_sim->Parameters().cornerMass *= cvol / 10e-15;


    _sim->InitVerticesAndElements(data->directors.get(), data->voxels, celldims);


    // Initialize the film graphics
    _film = Geometry{};

    // Function to access data
    std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* data, std::size_t i) {
        Magnum::Vector3* vdata = (Magnum::Vector3*)data;
        return vdata[i];
    };

    Vertex* nodes = _sim->getNodes();
    _numSurfaceNodes = 0;

    Vertex::Role role = Vertex::Role::Surface;

    for (int i = 0; i < _sim->getNumNodes(); i++) {
        if (nodes[i].role & role)
            ++_numSurfaceNodes;
    }

    _surfaceNodes.resize(_numSurfaceNodes);
    std::vector<Magnum::Vector3> tempNodes;
    tempNodes.resize(_numSurfaceNodes);

    // Create list of surface only nodes
    unsigned int count = 0;

    for (int i = 0; i < _sim->getNumNodes(); i++) {
        if (nodes[i].role & role) {

            float theta = acos(nodes[i].director[2]);
            float phi = atan2(nodes[i].director[1], nodes[i].director[0]);

            auto color = LC::Imaging::Colors::RungeSphere(theta, phi);

            _surfaceNodes[count] = i;
            // Convert to um
            tempNodes[count] = 1e6f*Vector3(nodes[i].position[0], nodes[i].position[1], nodes[i].position[2]);
            
            ++count;
        }
    }

    _film->Init(&tempNodes[0], access, _numSurfaceNodes);

    // Populate the z profile
    _zprofile = Zprofile{};
    
    updateFilmGraphics();
            
    updateBoxes();

    LC_INFO("Created client application!");
}

Sandbox::~Sandbox() {
    LC_INFO("Destroying client application.");
}

/*
    Main simulation loop + gui
*/
void Sandbox::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    _imgui.newFrame();
    
    mainWindow();

    /* Update application cursor */
    _imgui.updateApplicationCursor(*this);

    updateCamera();
}

Dataset* Sandbox::dataset() {
    return (Dataset*)_solver->GetDataPtr();
}


void Sandbox::updateCamera() {

    /* Update camera */
    if (_cameraType == CameraType::ArcBall) {
        _arcballCamera->updateTransformation();

        polyRenderer();
        updateBoxes();
        drawBoxes();

        if (_widget.updateFilm)
            updateFilm();

        if (_widget.drawFilm)
            drawFilm();

        if (_widget.drawZProfile)
            drawProfile();

        {
            guiRenderer();
            _imgui.drawFrame();
        }
        swapBuffers();

        redraw();
    }

}

void Sandbox::keyPressEvent(KeyEvent& event) {
    // Check if Ctrl + S or Ctrl + O is pressed
    if ((event.key() == KeyEvent::Key::S) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { save(); }
    else if ((event.key() == KeyEvent::Key::O) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) {
        if (open()) {
            // Change the state
        }
    }
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::updateBoxes(bool first) {

    _boxes.Clear();

    auto dims = _sim->Parameters().cell_dims;
    Vector3 s;
    
    s[0] = dims[0] * 0.5e6f;
    s[1] = dims[1] * 0.5e6f;
    s[2] = dims[2] * 0.5e6f;


    // Draw a 1x1x1 cube to the screen for reference
    _boxes.Append(_arcballCamera->viewMatrix() *
         Matrix4::scaling(s), 0x00ffff_rgbf);
}

void Sandbox::updateFilm() {

    for (int i = 0; i < _widget.updateCycle; i++)
        _sim->updateNodes();

    if (_widget.updateGraphics) {
        updateFilmGraphics();
    }
    // Film has been updated
    if (!_widget.continuousUpdate)
        _widget.updateFilm = false;
}


void Sandbox::updateFilmGraphics() {

    Dataset* data = dataset();
    Float polyRadius = _film->scale / (Float)pow(_film->numObjects, 1.0f / 3.0f);

    Matrix4 rotation;
    Float hPi = M_PI / 2.0;

    auto nodes = _sim->getNodes();
    auto voxels = _sim->Parameters().voxels;

    for (int i = 0; i < _numSurfaceNodes; i++) {

        float theta = acos(nodes[_surfaceNodes[i]].director[2]);
        float phi = atan2(nodes[_surfaceNodes[i]].director[1], nodes[_surfaceNodes[i]].director[0]);

        _film->polyPositions[i] = 1e6 * Vector3(nodes[_surfaceNodes[i]].position[0],
            nodes[_surfaceNodes[i]].position[1],
            nodes[_surfaceNodes[i]].position[2]);

        rotation = Matrix4::rotationZ(Rad{ -hPi + phi }) * Matrix4::rotationX(Rad{ hPi - theta });
        _film->polyInstanceData[i].transformationMatrix =
            Matrix4::translation(_film->polyPositions[i]) * Matrix4::scaling(Vector3{ polyRadius }) * rotation;
        _film->polyInstanceData[i].normalMatrix =
            _film->polyInstanceData[i].transformationMatrix.normalMatrix();
    }

    // Z profile
    std::unique_ptr<LC::scalar[]> rz(new LC::scalar[data->voxels[0] * data->voxels[1]]);
    for (int i = 0; i < voxels[0]; i++)
        for (int j = 0; j < voxels[1]; j++) {
            unsigned int idx = i + voxels[0] * j;
            rz[idx] = (_sim->getNodes() + (voxels[2] - 1) * voxels[0] * voxels[1] + idx)->position(2) * 1e6;
        }

    // Convert to um
    auto cdims = _sim->Parameters().cell_dims;
    for (auto &c : cdims)
        c *= 1e6;
    
    // Position rz has units um like cell_dims
    _zprofile->GenerateProfile(rz.get(), voxels, cdims);
}

void Sandbox::drawFilm() {
    _film->Draw(_arcballCamera, _projectionMatrix);
}

void Sandbox::drawProfile() {
    if (_zprofile)
        _zprofile->Draw(_arcballCamera, _projectionMatrix);
}

void Sandbox::drawBoxes() {
    _boxes.DrawFrame(_projectionMatrix);
}

void Sandbox::mainWindow() {

    if (!_widget.showSettings) return;

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::Begin("film activation", 0, window_flags);

    // Simulation details
    {
        std::string pitch = "pitch = " + std::to_string(_widget.pitch) + " um";
        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, pitch.c_str());
        std::string dims = "Cell dimensions = ("
            + std::to_string(_sim->Parameters().cell_dims[0] * 1e6)
            + " um X "
            + std::to_string(_sim->Parameters().cell_dims[1] * 1e6)
            + " um X "
            + std::to_string(_sim->Parameters().cell_dims[2] * 1e6)
            + " um)";
        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, dims.c_str());
    }

    ImGui::SliderFloat("Node scale", &_film->scale, 0.125f, 10.0f);
    ImGui::SliderInt("Chain extensibility", &_sim->Parameters().I_m, 4, 30);

    float kappa = _sim->Parameters().kappa * 1e-3;
    ImGui::SliderFloat("Bulk modulus (kPa)", &kappa, 0.f, 1e5f);
    _sim->Parameters().kappa = kappa * 1e3;

    float mu = _sim->Parameters().mu * 1e-3;
    ImGui::SliderFloat("Shear modulus (kPa)", &mu, 0.f, 800.f);
    _sim->Parameters().mu = mu * 1e3;

    float UdS = _sim->Parameters().UdS *1e-3;
    ImGui::SliderFloat("UdS (kPa)", &UdS, -800.f, 800.f);
    _sim->Parameters().UdS = UdS * 1e3;

    float dT = _sim->Parameters().dT * 1e9;
    ImGui::InputFloat("Time step (ns)", &dT);
    _sim->Parameters().dT = dT * 1e-9;

    ImGui::InputInt("Update iterations", &_widget.updateCycle);
    

    if (ImGui::Button("Update film"))
        _widget.updateFilm = true;

    ImGui::SameLine();

    ImGui::Checkbox("Continuous update", &_widget.continuousUpdate);

    ImGui::TextColored({ 1.f, 1.f, 0.f, 1.f }, "Draw");
    ImGui::Checkbox("Draw film", &_widget.drawFilm);
    ImGui::SameLine();
    ImGui::Checkbox("Draw z-profile", &_widget.drawZProfile);
    ImGui::Checkbox("Update graphics", &_widget.updateGraphics);

    // Show the film
    plotFilmData();


    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
        1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
    ImGui::End();
}

void Sandbox::plotFilmData() {
    Dataset* data = dataset();
    if (ImGui::RadioButton("x-line", &_widget.radioZProfileAxis, 0))
        _widget.radioZProfileIndex = -1;
    ImGui::SameLine();
    if (ImGui::RadioButton("y-line", &_widget.radioZProfileAxis, 1))
        _widget.radioZProfileIndex = -1;
    ImGui::SameLine();
    ImGui::PushItemWidth(150.0f);
    if (_widget.radioZProfileIndex == -1)
        _widget.radioZProfileIndex = (data->voxels[_widget.radioZProfileAxis] - 1) / 2;
    ImGui::SliderInt("Index", &_widget.radioZProfileIndex, 0, data->voxels[_widget.radioZProfileAxis] - 1);
    ImGui::PopItemWidth();

    if (_zprofile) {
        // Plot of the z-profile along x slice

        int vox[] = { _zprofile->graph->vox_x, _zprofile->graph->vox_x };
        int complementAxis = (_widget.radioZProfileAxis + 1) % 2;

        int vox_x = vox[_widget.radioZProfileAxis];
        int vox_y = vox[complementAxis];

        std::unique_ptr<float[]> ax(new float[vox_x]);
        std::unique_ptr<float[]> line_scan(new float[vox_y]);

        // Convert to um
        auto cdims = _sim->Parameters().cell_dims;
        for (auto& c : cdims)
            c *= 1e6;
        auto voxels = _sim->Parameters().voxels;



        float dx = cdims[_widget.radioZProfileAxis] / (voxels[_widget.radioZProfileAxis] - 1);
        float dy = cdims[complementAxis] / (voxels[complementAxis] - 1);

        for (int i = 0; i < vox_x; i++) {
            if (_widget.radioZProfileAxis == 0) {
                ax[i] = _zprofile->graph->data[i + _widget.radioZProfileIndex * vox[0]].x(); // um
                line_scan[i] = _zprofile->graph->data[i + _widget.radioZProfileIndex * vox[0]].z() * 1e3; // nm
            }
            else if (_widget.radioZProfileAxis == 1) {
                ax[i] = _zprofile->graph->data[_widget.radioZProfileIndex + i * vox[0]].y(); // um
                line_scan[i] = _zprofile->graph->data[_widget.radioZProfileIndex + i * vox[0]].z() * 1e3; // nm
            }
        }

        std::string letter;

        if (_widget.radioZProfileAxis == 0) // x-axis
            letter = 'x';
        else
            letter = 'y';

        std::string plt_title = "Z-profile plot (" + letter + "-axis)";

        auto fautofit = ImPlotAxisFlags_::ImPlotAxisFlags_AutoFit;

        if (ImPlot::BeginPlot(plt_title.c_str(), "L (um)", "h (nm)", ImVec2(-1, 0), 0, fautofit, fautofit)) {
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine(("Z(" + letter + ")").c_str(), ax.get(), line_scan.get(), vox_x);
            ImPlot::EndPlot();
        }

        // Draw the plane corresponding to the cross section

        _zprofile->box.dims[complementAxis] = 0;
        _zprofile->box.dims[_widget.radioZProfileAxis] = data->cell_dims[_widget.radioZProfileAxis];
        _zprofile->box.dims[2] = data->cell_dims[2];

        _zprofile->box.translation[complementAxis] = dy * (_widget.radioZProfileIndex - data->voxels[complementAxis] + 1) + 0.5f * data->cell_dims[complementAxis];
        _zprofile->box.translation[_widget.radioZProfileAxis] = 0.0f;
        _zprofile->box.translation[2] = 0.0f;

    }
}

LC::Application* LC::createApplication(int argc, char** argv) {

    return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}