#include <lclab2.h>
#include "Widget.h"

using namespace Magnum;
using namespace Math::Literals;
using AppSolver = LC::FrankOseen::ElasticOnly::RBFFDSolver;
using Dataset = AppSolver::Dataset;

struct BoxInstanceData {
    Matrix4 transformationMatrix;
    Color3 color;
};

class Sandbox : public LC::Application
{
public:
    explicit Sandbox(const Arguments& arguments);
    ~Sandbox();

    void keyPressEvent(KeyEvent& event) override;

    void initSpheres();


private:
    virtual void drawEvent() override;

    template <typename T>
    void dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map);

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader;
    SceneGraph::DrawableGroup3D _transparentDrawables;

    std::unique_ptr<LC::SphereArray> _sarray;

    GL::Mesh _boxMesh{ NoCreate };
    GL::Buffer _boxInstanceBuffer{ NoCreate };
    Shaders::FlatGL3D _boxShader{ NoCreate };
    Containers::Array<BoxInstanceData> _boxInstanceData;

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

    setupCamera(0.1f, CameraType::ArcBall);

    _solver = std::make_unique<AppSolver>();

    /* Setup data */
    Dataset* data = (Dataset*)_solver->GetDataPtr();
    (*data).ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        .Cell(1.0, 1.0, 1.0)
        .Boundaries(0, 0, 0)
        .Neighbors(16);

    // Things to add:
    // - Exclusion radius function presets
    // - Director configuration presets

    _solver->Init();

    //LC::Math::Weight<LC::scalar>* weight = data->derivative.GetWeight(LC::Math::WeightTag::x);

    //for (int i = 0; i < data->subnodes; i++) {

    //    std::string line;

    //    for (int k = 0; k < data->knn; k++) {
    //        line += std::to_string(weight->data[data->knn * i + k]) + " ";
     //   }
    //    LC_INFO("{0}", line.c_str());
    //}

    LC_INFO("Number of nodes = {0}", (*data).nodes);

    /* Init visuals */

    initSpheres();

    _boxShader = Shaders::FlatGL3D{
                    Shaders::FlatGL3D::Flag::VertexColor |
                    Shaders::FlatGL3D::Flag::InstancedTransformation };
    _boxInstanceBuffer = GL::Buffer{};
    _boxMesh = MeshTools::compile(Primitives::cubeWireframe());
    _boxMesh.addVertexBufferInstanced(_boxInstanceBuffer, 1, 0,
        Shaders::FlatGL3D::TransformationMatrix{},
        Shaders::FlatGL3D::Color3{});


    arrayAppend(_boxInstanceData, InPlaceInit,
        _arcballCamera->viewMatrix()*
        Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0], (float)data->cell_dims[1], (float)data->cell_dims[2] }), 0x00ffff_rgbf);

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
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    _imgui.newFrame();

    //Dataset* data = (Dataset*)(_solver->GetDataPtr());
    {
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;

        if (_widget.showSettings) {

            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::Begin("lclab2", 0, window_flags);

            {
                std::function<void()> menuFunc = [this]() {
                    initSpheres();
                };

                bool loaded = false;

                saveMenu(loaded, menuFunc);
            }



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
    if (_cameraType == CameraType::ArcBall) {
        _arcballCamera->updateTransformation();

        _boxInstanceData[0].transformationMatrix = _arcballCamera->viewMatrix() *
            Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0], (float)data->cell_dims[1], (float)data->cell_dims[2] });

        _boxInstanceBuffer.setData(_boxInstanceData, GL::BufferUsage::DynamicDraw);
        _boxMesh.setInstanceCount(_boxInstanceData.size());
        _boxShader.setTransformationProjectionMatrix(_projectionMatrix)
            .draw(_boxMesh);
    
        polyRenderer();
        _sarray->Draw(_arcballCamera, _projectionMatrix);

    }


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

void Sandbox::keyPressEvent(KeyEvent& event) {
    // Check if Ctrl + S or Ctrl + O is pressed
    if ((event.key() == KeyEvent::Key::S) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { save(); }
    else if ((event.key() == KeyEvent::Key::O) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { 
        if (open()) {
            initSpheres();
        }
    }
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::initSpheres() {
    // Spheres
    Dataset* data = (Dataset*)_solver->GetDataPtr();
    std::size_t Nspheres = (*data).nodes;
    std::function<Vector3(void*, std::size_t)> Identity = [Nspheres](void* data, std::size_t i) {
        LC::scalar* scalar_data = (LC::scalar*)data;
        return Vector3{ (float)scalar_data[i], (float)scalar_data[i + Nspheres],  (float)scalar_data[i + 2 * Nspheres] };
    };

    _sarray = std::unique_ptr <LC::SphereArray>(new LC::SphereArray);
    _sarray->Init((void*)(*data).position.get(), Identity, (*data).nodes, 2);

    for (int i = 0; i < data->subnodes; i++) {
        _sarray->sphereInstanceData[data->active_nodes[i]].color = Color3::red();
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