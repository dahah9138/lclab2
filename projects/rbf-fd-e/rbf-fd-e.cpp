#include <lclab2.h>
#include "Widget.h"

using namespace Magnum;
using namespace Math::Literals;
using AppSolver = LC::FrankOseen::ElasticOnly::RBFFDSolver;
using Dataset = AppSolver::Dataset;
using Geometry = LC::EllipsoidArray;

struct Plane {
    LC::Math::Interpolant<LC::scalar> interpolant;
    std::unique_ptr<LC::scalar[]> nodes;
    std::unique_ptr<std::size_t[]> neighbors;
    std::unique_ptr<Geometry> geometry;
    std::size_t numNodes;
    std::size_t knn = 4, xDensity = 20, yDensity = 20;
};


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

    void initGeometry();
    void updateGeometry();
    void generateInterpolant();
    void updateInterpolant();
    void updatePlane();


private:
    virtual void drawEvent() override;

    template <typename T>
    void dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map);

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader;
    SceneGraph::DrawableGroup3D _transparentDrawables;

    std::unique_ptr<Geometry> _garray;


    Plane _plane;

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
        .Boundaries(1, 1, 0)
        .Neighbors(16)
        .DirectorConfiguration(Dataset::Planar(1, 1.0));

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

    initGeometry();
    generateInterpolant();

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

    
    updatePlane();

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
                    initGeometry();
                    generateInterpolant();
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

            if (ImGui::Button("Update visuals")) {
                updateGeometry();
                updatePlane();
            }

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

        //_garray->Draw(_arcballCamera, _projectionMatrix);
        _plane.geometry->Draw(_arcballCamera, _projectionMatrix);

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
            initGeometry();
        }
    }
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::initGeometry() {

    Dataset* data = (Dataset*)_solver->GetDataPtr();
    std::size_t Nspheres = (*data).nodes;
    std::function<Vector3(void*, std::size_t)> Identity = [Nspheres](void* data, std::size_t i) {
        LC::scalar* scalar_data = (LC::scalar*)data;
        return Vector3{ (float)scalar_data[i], (float)scalar_data[i + Nspheres],  (float)scalar_data[i + 2 * Nspheres] };
    };

    _garray = std::unique_ptr <Geometry>(new Geometry);
    _garray->Init((void*)(*data).position.get(), Identity, (*data).nodes, 2);

    for (int i = 0; i < data->subnodes; i++) {
        _garray->polyInstanceData[data->active_nodes[i]].color = Color3::red();
    }

    updateGeometry();
}

void Sandbox::updateGeometry() {
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    Float polyRadius = _garray->scale / (Float)pow(_garray->numObjects, 1.0f / 3.0f);

    Matrix4 rotation;

    Float theta, phi;
    Float hPi = M_PI / 2.0;

    for (int i = 0; i < data->nodes; i++) {

        theta = acos(data->directors[i + 2 * data->nodes]);
        phi = atan2(data->directors[i + data->nodes], data->directors[i]);

        rotation = Matrix4::rotationZ(Rad{ -hPi + phi }) * Matrix4::rotationX(Rad{ hPi - theta });

        _garray->polyInstanceData[i].transformationMatrix = Matrix4::translation(_garray->polyPositions[i]) * Matrix4::scaling(Vector3{ polyRadius }) * rotation;
    }
}

void Sandbox::updatePlane() {

    updateInterpolant();
    
    Dataset* data = (Dataset*)_solver->GetDataPtr();
    std::size_t iX = _plane.xDensity * data->cell_dims[0];
    std::size_t iY = _plane.yDensity * data->cell_dims[1];

    LC::Math::Metric<LC::scalar> metric;
    metric.Bcs = data->bc;
    metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);

    Float theta, phi;
    Matrix4 rotation;
    Float hPi = M_PI / 2.0;

    Float polyRadius = _garray->scale / (Float)pow(_garray->numObjects, 1.0f / 3.0f);

    // Evaluate the interpolant for all points in plane
    for (int i = 0; i < iX; i++) {
        for (int j = 0; j < iX; j++) {

            std::size_t idx = j * iX + i;

            // Evaluate
            auto nn = _plane.interpolant.Evaluate(idx, data->position.get(), _plane.nodes.get(), _plane.neighbors.get(),
                *(data->RBF.get()), metric, _plane.numNodes, data->nodes, _plane.knn);

            // Normalize nn
            LC::scalar nmag = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
            nn[0] /= nmag;
            nn[1] /= nmag;
            nn[2] /= nmag;


            theta = acos(nn[2]);
            phi = atan2(nn[1], nn[0]);

            // Use runge sphere color
            _plane.geometry->polyInstanceData[idx].color = LC::Imaging::Colors::RungeSphere(theta, phi);

            // Set orientation
            rotation = Matrix4::rotationZ(Rad{ -hPi + phi }) * Matrix4::rotationX(Rad{ hPi - theta });
            _plane.geometry->polyInstanceData[idx].transformationMatrix = Matrix4::translation(_plane.geometry->polyPositions[idx]) * Matrix4::scaling(Vector3{ polyRadius }) * rotation;
        }
    }


}

void Sandbox::generateInterpolant() {
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    // 1. Create the node positions
    std::size_t iX = _plane.xDensity * data->cell_dims[0];
    std::size_t iY = _plane.yDensity * data->cell_dims[1];

    _plane.numNodes = iX * iY;
    _plane.nodes = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * _plane.numNodes]);
    _plane.neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[_plane.knn * _plane.numNodes]);
    


    for (int i = 0; i < iX; i++) {
        for (int j = 0; j < iY; j++) {
            std::size_t idx = j * iX + i;
            LC::scalar xx = (LC::scalar)i / LC::scalar(iX - 1) - 0.5;
            LC::scalar yy = (LC::scalar)j / LC::scalar(iY - 1) - 0.5;

            _plane.nodes[idx] = xx * data->cell_dims[0];
            _plane.nodes[idx + _plane.numNodes] = yy * data->cell_dims[1];
            _plane.nodes[idx + 2 * _plane.numNodes] = 0.0;
        }
    }

    // Initialize the plane geometry
    {
        std::size_t Nspheres = _plane.numNodes;
        std::function<Vector3(void*, std::size_t)> Identity = [Nspheres](void* data, std::size_t i) {
            LC::scalar* scalar_data = (LC::scalar*)data;
            return Vector3{ (float)scalar_data[i], (float)scalar_data[i + Nspheres],  (float)scalar_data[i + 2 * Nspheres] };
        };

        _plane.geometry = std::unique_ptr <Geometry>(new Geometry);
        _plane.geometry->Init((void*)_plane.nodes.get(), Identity, _plane.numNodes, 2);
    }

    LC::Math::Metric<LC::scalar> metric;
    metric.Bcs = data->bc;
    metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);

    // 2. Compute neighbors relative to actual node positions
    LC::Algorithm::knn_c(data->position.get(), data->nodes, _plane.nodes.get(), _plane.numNodes, metric, _plane.knn, (LC::scalar*)0, _plane.neighbors.get());

    // 3. Compute interpolant
    _plane.interpolant.ComputeFactorization(data->position.get(), _plane.neighbors.get(), *(data->RBF.get()), metric, _plane.numNodes, data->nodes, _plane.knn);
}

void Sandbox::updateInterpolant() {
    Dataset* data = (Dataset*)_solver->GetDataPtr();
    _plane.interpolant.ComputeWeights(data->directors.get(), data->neighbors.get(), _plane.numNodes, data->nodes, _plane.knn);
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