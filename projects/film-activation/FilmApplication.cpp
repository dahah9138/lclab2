#include "DropDownMenu.h"
#include "Widget.h"
#include "Boxes.h"
#include "Simulation.h"

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
    void mainWindow();
    void updateCamera();
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

    _sim = std::make_unique<Simulation>();
    
    // Default simulation initialization
    _sim->InitVerticesAndElements();

    // Initialize the film graphics
    _film = Geometry{};

    // Function to access data
    std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* data, std::size_t i) {
        Magnum::Vector3* vdata = (Magnum::Vector3*)data;
        return vdata[i];
    };

    Vertex* nodes = _sim->getNodes();
    _numSurfaceNodes = 0;

    Vertex::Role role = Vertex::Role::Bulk;

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
            tempNodes[count] = Vector3(nodes[i].position[0], nodes[i].position[1], nodes[i].position[2]);
            
            ++count;
        }
    }

    _film->Init(&tempNodes[0], access, _numSurfaceNodes);
    
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

        drawFilm();

        // Sort objects to draw in correct order

        //sortObjects(_transparentDrawables);
        //_camera->draw(_transparentDrawables);
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

    // Draw a 1x1x1 cube to the screen for reference
    _boxes.Append(_arcballCamera->viewMatrix() *
        Matrix4::translation(Vector3{ 0.f, 0.f, -0.8f }) * Matrix4::scaling(Vector3{ 1.f,1.f,.2f }), 0x00ffff_rgbf);
}

void Sandbox::updateFilm() {

    for (int i = 0; i < _widget.updateCycle; i++)
        _sim->updateNodes();

    updateFilmGraphics();
    // Film has been updated
    _widget.updateFilm = false;
}

void Sandbox::updateFilmGraphics() {
    // Update the film (?)
    Float polyRadius = _film->scale / (Float)pow(_film->numObjects, 1.0f / 3.0f);

    Matrix4 rotation;
    Float hPi = M_PI / 2.0;

    auto nodes = _sim->getNodes();

    for (int i = 0; i < _numSurfaceNodes; i++) {

        float theta = acos(nodes[_surfaceNodes[i]].director[2]);
        float phi = atan2(nodes[_surfaceNodes[i]].director[1], nodes[_surfaceNodes[i]].director[0]);

        _film->polyPositions[i] = Vector3(nodes[_surfaceNodes[i]].position[0],
            nodes[_surfaceNodes[i]].position[1],
            nodes[_surfaceNodes[i]].position[2]);

        rotation = Matrix4::rotationZ(Rad{ -hPi + phi }) * Matrix4::rotationX(Rad{ hPi - theta });
        _film->polyInstanceData[i].transformationMatrix =
            Matrix4::translation(_widget.positionScale * _film->polyPositions[i]) * Matrix4::scaling(Vector3{ polyRadius }) * rotation;
        _film->polyInstanceData[i].normalMatrix =
            _film->polyInstanceData[i].transformationMatrix.normalMatrix();
    }
}

void Sandbox::drawFilm() {
    _film->Draw(_arcballCamera, _projectionMatrix);
}

void Sandbox::drawBoxes() {
    _boxes.DrawFrame(_projectionMatrix);
}

void Sandbox::mainWindow() {

    if (!_widget.showSettings) return;

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::Begin("lclab2", 0, window_flags);

    ImGui::SliderFloat("Node scale", &_film->scale, 0.125f, 10.0f);
    ImGui::SliderFloat("Position scale", &_widget.positionScale, 0.005f, 10.f);
    ImGui::InputInt("Update iterations", &_widget.updateCycle);
    //if (ImGui::Button("Update film"))
        _widget.updateFilm = true;

    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
        1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
    ImGui::End();
}


LC::Application* LC::createApplication(int argc, char** argv) {

    return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}