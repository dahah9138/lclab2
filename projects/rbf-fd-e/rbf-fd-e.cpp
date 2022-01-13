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
    std::size_t knn = 35, xDensity = 10, yDensity = 10;
    // Track usability
    bool usable = false;
    bool generated = false;
};

struct Region {
    std::unique_ptr<Geometry> geometry;
    std::vector<std::pair<std::size_t, std::array<float,3>>> nodes;
    std::array<float, 3> interval = { 1.0f, 1.0f, 0.1f };
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
    void updateBoxes(bool first = false);
    void drawBoxes();

private:
    virtual void drawEvent() override;

    template <typename T>
    void dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map);

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader;
    SceneGraph::DrawableGroup3D _transparentDrawables;
    
    Region _region;

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

    LC::scalar dop = 8.0;
    int Q = 2;
    LC::scalar qRatio = dop / (2. * Q + 1.);

    LC::scalar active_rad = dop / (2. * qRatio);
    LC::scalar padding = dop * 0.1;
    LC::scalar r_edge = dop / 2.0;

    (*data).ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        .Cell(dop, dop, dop)
        .Boundaries(0, 0, 0)
        .Neighbors(80)
        .IsActiveConfig(LC::Math::ActiveSphere(active_rad))
        .ExclusionRadius(LC::Math::LinearSphere(0.2, 0.25, r_edge, active_rad + padding))
        //.ExclusionRadius(Dataset::UniformRadius(0.15))
        .DirectorConfiguration(LC::Math::Heliknoton(Q, { dop, dop, dop }, 1. / qRatio));
        //.DirectorConfiguration(Dataset::Planar(2 * dop, dop));
        //.DirectorConfiguration(Dataset::Uniform());
    _solver->Init();

    LC_INFO("Number of nodes = {0}", (*data).nodes);

    /* Init visuals */

    if (_widget.interpolant) updatePlane();
    else initGeometry();

    _boxShader = Shaders::FlatGL3D{
                    Shaders::FlatGL3D::Flag::VertexColor |
                    Shaders::FlatGL3D::Flag::InstancedTransformation };
    _boxInstanceBuffer = GL::Buffer{};
    _boxMesh = MeshTools::compile(Primitives::cubeWireframe());
    _boxMesh.addVertexBufferInstanced(_boxInstanceBuffer, 1, 0,
        Shaders::FlatGL3D::TransformationMatrix{},
        Shaders::FlatGL3D::Color3{});

    updateBoxes(true);


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
                    if (_widget.interpolant) {
                        generateInterpolant();
                        updatePlane();
                    }
                    else initGeometry();
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
                float scale = 50.0f * _widget.nodeScale;
                ImGui::SliderFloat("Node scale", &scale, 0.01f, 4.0f);
                _widget.nodeScale = scale / 50.0f;
            }

            {
                bool prev = _widget.interpolant;
                ImGui::Checkbox("Interpolant", &_widget.interpolant);

                if (!prev && _widget.interpolant) {
                    _plane.generated = false;
                    _plane.usable = false;
                    _widget.updateImage = true;
                }
                else if (prev && !_widget.interpolant) {
                    _widget.updateImage = true;
                }
            }

            if (!_widget.interpolant) {
                ImGui::SliderFloat("X ratio", &_region.interval[0], 0.0f, 1.0f);
                ImGui::SliderFloat("Y ratio", &_region.interval[1], 0.0f, 1.0f);
                ImGui::SliderFloat("Z ratio", &_region.interval[2], 0.0f, 1.0f);
            }
            else {

            }

            // Pressed the relax button
            ImGui::InputInt("Iterations/cycle", &_widget.cycle);
            _widget.relax = ImGui::Button("Relax");
            if (_relaxFuture.second) {
                ImGui::SameLine();
                ImGui::Text("Relaxing...");
            }

            if (ImGui::Button("Update visuals")) {
                initGeometry();
            }

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));

            ImGui::End();

        }


        if (_widget.relax) {

            const bool GPU = true;

            // Try to relax asynchronously...
            _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, GPU);
            _relaxFuture.second = true;
            _widget.updateImage = true;
        }

        if (_widget.updateImage || _relaxFuture.second) {

            bool ready = checkRelax();

            if (ready) {

                if (_widget.interpolant) updatePlane();
                else initGeometry();

                _widget.updateImage = false;
            }

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

    
        polyRenderer();

        updateBoxes();
        drawBoxes();

        if (_widget.interpolant) _plane.geometry->Draw(_arcballCamera, _projectionMatrix);
        else  _region.geometry->Draw(_arcballCamera, _projectionMatrix);
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
            if (_widget.interpolant) {
                generateInterpolant(); 
                updatePlane();
            }
            else initGeometry();
        }
    }
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::initGeometry() {

    Dataset* data = (Dataset*)_solver->GetDataPtr();

    {
        std::vector<std::pair<std::size_t, std::array<float, 3>>> temp{};
        _region.nodes.swap(temp);
        _region.nodes.reserve(data->nodes);
    }
    

    // Add points within bounds
    for (std::size_t i = 0; i < data->nodes; i++) {
        for (int d = 0; d < 3; d++) {

            float bound = _region.interval[d]  * data->cell_dims[d] * 0.5;

            float pd = abs(data->position[i + data->nodes * d]);

            if (pd > bound)
                break;
            
            // valid node
            if (d == 2) {
                std::array<float, 3> pos{ data->position[i], data->position[i + data->nodes], data->position[i + 2 * data->nodes] };
                _region.nodes.emplace_back(i, pos);
            }
        }
    }

    std::size_t Nshapes = _region.nodes.size();

    std::function<Vector3(void*, std::size_t)> Identity = [Nshapes](void* data, std::size_t i) {
        std::pair<std::size_t, std::array<float, 3>>* p_data = (std::pair<std::size_t, std::array<float, 3>>*)data;
        return Vector3{ p_data[i].second[0], p_data[i].second[1], p_data[i].second[2] };
    };

    // Initialize the geometry
    _region.geometry = std::unique_ptr <Geometry>(new Geometry);
    _region.geometry->Init((void*)&_region.nodes[0], Identity, _region.nodes.size(), 2);

    updateGeometry();
}

void Sandbox::updateGeometry() {
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    Matrix4 rotation;

    Float theta, phi;
    Float hPi = M_PI / 2.0;

    for (std::size_t i = 0; i < _region.geometry->numObjects; i++) {

        std::size_t idx = _region.nodes[i].first;

        theta = acos(data->directors[idx + 2 * data->nodes]);
        phi = atan2(data->directors[idx + data->nodes], data->directors[idx]);

        rotation = Matrix4::rotationZ(Rad{ -hPi + phi }) * Matrix4::rotationX(Rad{ hPi - theta });

        _region.geometry->polyInstanceData[i].transformationMatrix = Matrix4::translation(_region.geometry->polyPositions[i]) * Matrix4::scaling(Vector3{ _widget.nodeScale }) * rotation;
        _region.geometry->polyInstanceData[i].color = LC::Imaging::Colors::RungeSphere(theta, phi);
    }
}

void Sandbox::updatePlane() {

    // Need to recompute weights
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

    // Evaluate the interpolant for all points in plane
    for (int i = 0; i < iX; i++) {
        for (int j = 0; j < iY; j++) {

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
            _plane.geometry->polyInstanceData[idx].transformationMatrix = Matrix4::translation(_plane.geometry->polyPositions[idx]) * Matrix4::scaling(Vector3{ _widget.nodeScale }) * rotation;
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
    
    // Generate the xy midplane
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

        _plane.geometry = std::unique_ptr<Geometry>(new Geometry);
        _plane.geometry->Init((void*)_plane.nodes.get(), Identity, _plane.numNodes, 2);
    }

    LC::Math::Metric<LC::scalar> metric;
    metric.Bcs = data->bc;
    metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);

    // 2. Compute neighbors relative to actual node positions
    LC::Algorithm::knn_c(data->position.get(), data->nodes, _plane.nodes.get(), _plane.numNodes, metric, _plane.knn, (LC::scalar*)0, _plane.neighbors.get());

    // 3. Compute interpolant
    _plane.interpolant.ComputeFactorization(data->position.get(), _plane.neighbors.get(), *(data->RBF.get()), metric, _plane.numNodes, data->nodes, _plane.knn);

    _plane.generated = true;
}

void Sandbox::updateInterpolant() {

    if (!_plane.generated) generateInterpolant();

    Dataset* data = (Dataset*)_solver->GetDataPtr();
    _plane.interpolant.ComputeWeights(data->directors.get(), _plane.neighbors.get(), _plane.numNodes, data->nodes, _plane.knn);
    _plane.usable = true;
}

void Sandbox::updateBoxes(bool first) {
    Dataset* data = (Dataset*)_solver->GetDataPtr();
    
    if (!first) arrayResize(_boxInstanceData, 0);

    arrayAppend(_boxInstanceData, InPlaceInit,
        _arcballCamera->viewMatrix() *
        Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0], (float)data->cell_dims[1], (float)data->cell_dims[2] }), 0x00ffff_rgbf);

    if (!_widget.interpolant) {
        arrayAppend(_boxInstanceData, InPlaceInit,
            _arcballCamera->viewMatrix() *
            Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0] * _region.interval[0],
                (float)data->cell_dims[1] * _region.interval[1],
                (float)data->cell_dims[2] * _region.interval[2] }), 0x00ffff_rgbf);
    }
}

void Sandbox::drawBoxes() {
    _boxInstanceBuffer.setData(_boxInstanceData, GL::BufferUsage::DynamicDraw);
    _boxMesh.setInstanceCount(_boxInstanceData.size());
    _boxShader.setTransformationProjectionMatrix(_projectionMatrix)
        .draw(_boxMesh);
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