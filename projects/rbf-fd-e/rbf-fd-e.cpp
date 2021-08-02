#include <lclab2.h>
#include "Widget.h"

using namespace Magnum;
using namespace Math::Literals;
using AppSolver = LC::FrankOseen::ElasticOnly::RBFFDSolver;
using Dataset = AppSolver::Dataset;

class Sandbox : public LC::Application
{
public:
    explicit Sandbox(const Arguments& arguments);
    ~Sandbox();

private:
    virtual void drawEvent() override;
    void keyPressEvent(KeyEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;

    template <typename T>
    void dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map);

    void save();
    bool open();
    void saveAs();

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader;
    SceneGraph::DrawableGroup3D _transparentDrawables;

    std::unique_ptr<LC::SphereArray> _sarray;

    Widget _widget;
    LC::Imaging::UniformGrid::POM _pomImager;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("FOFD Elastic Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable) } {
    _transparentShader = Shaders::VertexColorGL3D{};

    /* Setup the GUI */
    setupGUI();

    /* Setup window and parameters */
    enableDepthTest();
    enableFaceCulling();

    /* Loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);

    setupCamera(2.0f, CameraType::Group);

    _solver = std::make_unique<AppSolver>();

    /* Setup data */
    // TODO

    /* Init visuals */

    // TODO
    //_sarray.Init();
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

    //Dataset* data = (Dataset*)(_solver->GetDataPtr());
    {
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;

        if (_widget.showSettings) {

            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::Begin("lclab2", 0, window_flags);

            // Demo widgets for fast prototyping
            if (ImGui::Button("Demo ImGui"))
                _widget.showDemoWindow ^= true;

            ImGui::SameLine();

            static bool p_open = false;
            if (ImGui::Button("Demo ImPlot"))
                p_open ^= true;
            if (p_open)
                ImPlot::ShowDemoWindow(&p_open);

            if (ImGui::Button("LC Info"))
                _widget.showAnotherWindow ^= true;

            {
                ImGui::SliderFloat("Alpha", &_widget.alpha, 0.0f, 1.0f);
            }

            // Pressed the relax button
            _widget.relax = ImGui::Button("Relax");

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));

            ImGui::End();
        }
    }


    /* 3. Show the ImGui demo window. Most of the sample code is in
       ImGui::ShowDemoWindow() */
    if (_widget.showDemoWindow) {
        ImGui::SetNextWindowPos(ImVec2(100, 20), ImGuiCond_FirstUseEver);
        ImGui::ShowDemoWindow();
    }

    /* Update application cursor */
    _imgui.updateApplicationCursor(*this);

    /* Update camera */
    if (_cameraType == CameraType::ArcBall)
        _arcballCamera->updateTransformation();

    /* Reset state. Only needed if you want to draw something else with
        different state after. */

    polyRenderer();

    // Sort objects to draw in correct order
    sortObjects(_transparentDrawables);

    _camera->draw(_transparentDrawables);
    {
        guiRenderer();
        _imgui.drawFrame();
    }
    swapBuffers();

    redraw();
}

// Todo: inherit all these events through Application

void Sandbox::keyPressEvent(KeyEvent& event) {

    // Check if Ctrl + S or Ctrl + O is pressed


    if ((event.key() == KeyEvent::Key::S) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { /* ctrl + S action */}
    else if ((event.key() == KeyEvent::Key::O) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { /* ctrl + O action */ }
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::keyReleaseEvent(KeyEvent& event) {

    if ((event.key() == KeyEvent::Key::S)) {}
    else if ((event.key() == KeyEvent::Key::O)) {}
    if (_imgui.handleKeyReleaseEvent(event)) _ioUpdate = true;
}

void Sandbox::save() {

    bool ready = checkRelax();
    if (!ready) return;

    // Check if file exists (was loaded)

    std::string res;


    if (!_widget.loadedFromFile) {
        res = saveDialog();
    }
    else {
        res = _header.readFile;
    }

    if (!res.empty()) {
        AppSolver* solver = (AppSolver*)_solver.get();

        // Write ConfigureHeader in RBFFDSolver::Dataset
        //solver->ConfigureHeader(_header);

        _header.write(res);
        _header.writeBody();
        LC_INFO("Saved to file {0}", res.c_str());
        // If the file was saved, then it is now technically the same state
        // as 'loadedFromFile'
        _widget.loadedFromFile = true;
        _header.readFile = res;
    }

}

bool Sandbox::open() {

    bool ready = checkRelax();
    if (!ready) return false;
    bool updateImageFromLoad = false;

    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    auto res = openDialog();

    if (!res.empty()) {
        AppSolver* solver = (AppSolver*)_solver.get();

        _header.readFile = res[0];

        // Write ReadDataFromHeader in RBFFDSolver::Dataset
        //solver->ReadDataFromHeader(_header);

        updateImageFromLoad = true;

        // Tells program to save to readFile if "Save" is pressed
        _widget.loadedFromFile = true;
        LC_INFO("Loaded file {0}", res[0].c_str());

    }

    return updateImageFromLoad;
}

void Sandbox::saveAs() {

    bool ready = checkRelax();
    if (!ready) return;

    // Pass file to header to be written
    auto res = saveDialog();

    if (!res.empty()) {
        AppSolver* solver = (AppSolver*)_solver.get();
        
        // Write ConfigureHeader in RBFFDSolver::Dataset
        //solver->ConfigureHeader(_header);

        _header.write(res);
        _header.writeBody();
        LC_INFO("Saved to file {0}", res.c_str());
    }
}

template <typename T>
void Sandbox::dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map) {

    std::string currentItem = map[currentSelectable];

    if (ImGui::BeginCombo(menuName, currentItem.c_str())) {
        for (const auto& elem : map) {
            bool selected = (currentSelectable == elem.first);
            if (ImGui::Selectable(elem.second.c_str(), selected)) {
                currentSelectable = elem.first;
            }
            if (selected)
                ImGui::SetItemDefaultFocus();
        }

        ImGui::EndCombo();
    }
}

LC::Application* LC::createApplication(int argc, char** argv) {

    return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}