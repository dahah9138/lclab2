#include <lclab2.h>
#include "Widget.h"

#define TAYLOR_METHOD 0
#define HOPFION_CONFIG 0

using namespace Magnum;
using namespace Math::Literals;
using AppSolver = LC::FrankOseen::Electric::RBFFDSolver;
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

struct TaylorPlane {
    LC::Math::TaylorSeries<LC::scalar> series;
    std::unique_ptr<LC::scalar[]> nodes;
    std::unique_ptr<std::size_t[]> neighbors;
    std::unique_ptr<Geometry> geometry;
    std::size_t numNodes;
    std::size_t xDensity = 20, yDensity = 20;
    // Track usability
    bool usable = false;
};

struct Region {
    std::unique_ptr<Geometry> geometry;
    std::vector<std::pair<std::size_t, std::array<float, 3>>> nodes;
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
    TaylorPlane _tplane;

    GL::Mesh _boxMesh{ NoCreate };
    GL::Buffer _boxInstanceBuffer{ NoCreate };
    Shaders::FlatGL3D _boxShader{ NoCreate };
    Containers::Array<BoxInstanceData> _boxInstanceData;

    Widget _widget;
    LC::Imaging::UniformGrid::POM _pomImager;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("FO-RBF-FD Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable) } {
    _transparentShader = Shaders::VertexColorGL3D{};


    Utility::Arguments args;
    args.addOption('q', "topological-charge", "1")
        .setHelp("topological-charge", "topological charge", "Q")
        .addOption('k', "nearest-neighbors", "50")
        .setHelp("nearest-neighbors", "nearest neighbors", "K")
        .addOption('r', "exclusion-radius", "0.14285714285")
        .setHelp("exclusion-radius", "exclusion radius", "R")
        .addOption('v', "voltage", "2.0")
        .setHelp("voltage", "Voltage", "V")
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

    _solver = std::make_unique<AppSolver>();

    /* Setup data */
    Dataset* data = (Dataset*)_solver->GetDataPtr();

#if HOPFION_CONFIG

    int Q = 1;
    LC::scalar side = 2. * Q + 1.;
    LC::scalar voltage = 1.3;
    LC::scalar cellZ = 0.93;

    int knn = 80;
    LC::scalar rNode = 1. / 10.;
    data->npp = 30.;
    LC::scalar sigma = 1.;
    LC::scalar active_rad = side / 2;
    LC::scalar padding = 0;// pow(3. * knn / (4. * M_PI), 1. / 3.)* rNode;
    LC::scalar r_edge = side / 2.0;

    (*data)
        .ElectricConstants(LC::FrankOseen::LC_TYPE::_5CB)
        .VoltageConfiguration(LC::Math::VoltageZ(0.0, voltage, cellZ))
        .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        .Rate(-.50)
        .Cell(side, side, cellZ)
        .Boundaries(0, 0, 0)
        .Neighbors(knn)
        .ExclusionRadius(LC::Math::UniformRadius(rNode))
        .DirectorConfiguration(LC::Math::Hopfion(Q, { side, side, cellZ }, 1.0, 1.5));

    _region.interval[0] = 1.0;
    _region.interval[1] = 1.0;
#else
    /*
    Q = 1 : Stable (dop = 3, V = 2)
    Q = 2 : Stable (dop = 5, V = 2.7-2.8)
    Q = 3 : Stable (dop = 7, V = 3.13)
    */
    int Q = args.value<int>("topological-charge");
    LC::scalar voltage = args.value<LC::scalar>("voltage");
    LC::scalar dop = 2. * Q + 1.;

    int knn = args.value<int>("nearest-neighbors");
    LC::scalar rNode = args.value<LC::scalar>("exclusion-radius");
    data->npp = 30.;
    LC::scalar qRatio = dop / (2. * Q + 1.);
    LC::scalar sigma = 1.;
    LC::scalar active_rad = Q + 0.5;
    LC::scalar padding = 0;// pow(3. * knn / (4. * M_PI), 1. / 3.)* rNode;
    LC::scalar r_edge = dop / 2.0;

    (*data)
        .ElectricConstants(LC::FrankOseen::LC_TYPE::_5CB)
        .VoltageConfiguration(LC::Math::VoltageZ(0.0, voltage, dop))
        .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        .Rate(_widget.relaxRate)
        .Cell(dop, dop, dop)
        .Boundaries(0, 0, 0)
        .Neighbors(knn)
        //.IsActiveConfig(LC::Math::ActiveSphere(active_rad))
        //.ExclusionRadius(LC::Math::LinearSphere(rNode, sqrt(3.0) * rNode, sqrt(2.) * active_rad, active_rad))
        //.ExclusionRadius(LC::Math::ZLine(rNode, 0.25, 0.5, 1.5))
        .ExclusionRadius(LC::Math::UniformRadius(rNode))
        .DirectorConfiguration(LC::Math::Heliknoton(Q, { dop, dop, dop }, 1. / qRatio));
    //.DirectorConfiguration(Dataset::Planar(2 * dop, dop));
    //.DirectorConfiguration(Dataset::Uniform());

    _region.interval[0] = 2.0 * (active_rad +  padding) / dop;
    _region.interval[1] = 2.0 * (active_rad +  padding) / dop;
#endif

    _solver->Init();


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

#if TAYLOR_METHOD
                    _tplane.usable = false;
#else
                    _plane.generated = false;
                    _plane.usable = false;
#endif


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

                float intervalX = _region.interval[0];
                float intervalY = _region.interval[1];

                ImGui::SliderFloat("X ratio", &_region.interval[0], 0.0f, 1.0f);
                ImGui::SliderFloat("Y ratio", &_region.interval[1], 0.0f, 1.0f);

                if (intervalX != _region.interval[0] || intervalY != _region.interval[1]) {
                    // Need to regenerate the interpolant
                    _widget.regenerateInterpolant = true;
                }
            }

            // Pressed the relax button
            ImGui::PushItemWidth(100.0f);

            ImGui::InputInt("Iterations/cycle", &_widget.cycle);

            ImGui::SameLine();

            ImGui::InputFloat("Relax rate", &_widget.relaxRate);
            _widget.relax = ImGui::Button("Relax");


            ImGui::PopItemWidth();
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
            data->Rate(_widget.relaxRate);
            _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, GPU);

            _relaxFuture.second = true;
            _widget.updateImage = true;
        }

        if (_widget.updateImage || _relaxFuture.second) {

            bool ready = checkRelax();

            if (ready) {

                if (_widget.interpolant) {

                    if (_widget.regenerateInterpolant) {
                        generateInterpolant();
                        updateInterpolant();
                        _widget.regenerateInterpolant = false;
                    }

                    updatePlane();
                }
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

#if TAYLOR_METHOD
        if (_widget.interpolant) _tplane.geometry->Draw(_arcballCamera, _projectionMatrix);
#else
        if (_widget.interpolant) _plane.geometry->Draw(_arcballCamera, _projectionMatrix);
#endif
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

            float bound = _region.interval[d] * data->cell_dims[d] * 0.5;

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

    float cxRatio = data->cell_dims[0] * _region.interval[0];
    float cyRatio = data->cell_dims[1] * _region.interval[1];

#if TAYLOR_METHOD
    std::size_t iX = _tplane.xDensity * cxRatio;
    std::size_t iY = _tplane.yDensity * cyRatio;
#else
    std::size_t iX = _plane.xDensity * cxRatio;
    std::size_t iY = _plane.yDensity * cyRatio;
#endif

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
            //auto nn = _plane.interpolant.Evaluate(idx, data->position.get(), _plane.nodes.get(), _plane.neighbors.get(),
            //    *(data->RBF.get()), metric, _plane.numNodes, data->nodes, _plane.knn);
#if TAYLOR_METHOD
            auto nn = _tplane.series.Evaluate(idx, _tplane.nodes.get(), data->position.get(), data->directors.get(), _tplane.neighbors.get(), data->nodes);
#else
            auto nn = _plane.interpolant.Evaluate(idx, data->position.get(), _plane.nodes.get(), _plane.neighbors.get(),
                    *(data->RBF.get()), metric, _plane.numNodes, data->nodes, _plane.knn);
#endif
            // Normalize nn
            LC::scalar nmag = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
            nn[0] /= nmag;
            nn[1] /= nmag;
            nn[2] /= nmag;


            theta = acos(nn[2]);
            phi = atan2(nn[1], nn[0]);

            // Use runge sphere color
#if TAYLOR_METHOD
            _tplane.geometry->polyInstanceData[idx].color = LC::Imaging::Colors::RungeSphere(theta, phi);
#else
            _plane.geometry->polyInstanceData[idx].color = LC::Imaging::Colors::RungeSphere(theta, phi);
#endif
            // Set orientation
            rotation = Matrix4::rotationZ(Rad{ -hPi + phi }) * Matrix4::rotationX(Rad{ hPi - theta });

#if TAYLOR_METHOD
            _tplane.geometry->polyInstanceData[idx].transformationMatrix = Matrix4::translation(_tplane.geometry->polyPositions[idx]) * Matrix4::scaling(Vector3{ _widget.nodeScale }) * rotation;
#else
            _plane.geometry->polyInstanceData[idx].transformationMatrix = Matrix4::translation(_plane.geometry->polyPositions[idx]) * Matrix4::scaling(Vector3{ _widget.nodeScale }) * rotation;
#endif
        }
    }

}

void Sandbox::generateInterpolant() {
    Dataset* data = (Dataset*)_solver->GetDataPtr();


    // 1. Create the node positions

    float cxRatio = data->cell_dims[0] * _region.interval[0];
    float cyRatio = data->cell_dims[1] * _region.interval[1];

#if TAYLOR_METHOD
    std::size_t iX = _tplane.xDensity * cxRatio;
    std::size_t iY = _tplane.yDensity * cyRatio;
#else
    _plane.usable = false;
    std::size_t iX = _plane.xDensity * cxRatio;
    std::size_t iY = _plane.yDensity * cyRatio;
#endif

#if TAYLOR_METHOD
    _tplane.numNodes = iX * iY;
    _tplane.nodes = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * _tplane.numNodes]);
    _tplane.neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[_tplane.numNodes]);
#else
    _plane.numNodes = iX * iY;
    _plane.nodes = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * _plane.numNodes]);
    _plane.neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[_plane.knn * _plane.numNodes]);
#endif

    // Generate the xy midplane
    for (int i = 0; i < iX; i++) {
        for (int j = 0; j < iY; j++) {
            std::size_t idx = j * iX + i;
            LC::scalar xx = (LC::scalar)i / LC::scalar(iX - 1) - 0.5;
            LC::scalar yy = (LC::scalar)j / LC::scalar(iY - 1) - 0.5;

#if TAYLOR_METHOD
            _tplane.nodes[idx] = xx * cxRatio;
            _tplane.nodes[idx + _tplane.numNodes] = yy * cyRatio;
            _tplane.nodes[idx + 2 * _tplane.numNodes] = 0.0;
#else
            _plane.nodes[idx] = xx * cxRatio;
            _plane.nodes[idx + _plane.numNodes] = yy * cyRatio;
            _plane.nodes[idx + 2 * _plane.numNodes] = 0.0;
#endif
        }
    }

    // Initialize the plane geometry
    {
#if TAYLOR_METHOD
        std::size_t Nspheres = _tplane.numNodes;
#else
        std::size_t Nspheres = _plane.numNodes;
#endif
        std::function<Vector3(void*, std::size_t)> Identity = [Nspheres](void* data, std::size_t i) {
            LC::scalar* scalar_data = (LC::scalar*)data;
            return Vector3{ (float)scalar_data[i], (float)scalar_data[i + Nspheres],  (float)scalar_data[i + 2 * Nspheres] };
        };
#if TAYLOR_METHOD
        _tplane.geometry = std::unique_ptr<Geometry>(new Geometry);
        _tplane.geometry->Init((void*)_tplane.nodes.get(), Identity, _tplane.numNodes, 2);
#else
        _plane.geometry = std::unique_ptr<Geometry>(new Geometry);
        _plane.geometry->Init((void*)_plane.nodes.get(), Identity, _plane.numNodes, 2);
#endif
    }

    LC::Math::Metric<LC::scalar> metric;
    metric.Bcs = data->bc;
    metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);

    // 2. Get first nearest neighbor
#if TAYLOR_METHOD
    LC::Algorithm::knn_c(data->position.get(), data->nodes, _tplane.nodes.get(), _tplane.numNodes, metric, 1, (LC::scalar*)0, _tplane.neighbors.get());
#else
    LC::Algorithm::knn_c(data->position.get(), data->nodes, _plane.nodes.get(), _plane.numNodes, metric, _plane.knn, (LC::scalar*)0, _plane.neighbors.get());
#endif
    // 3. Compute interpolant
#if TAYLOR_METHOD
    LC::scalar* dx = data->derivative.GetWeight(LC::Math::WeightTag::x)->data;
    LC::scalar* dy = data->derivative.GetWeight(LC::Math::WeightTag::y)->data;
    LC::scalar* dz = data->derivative.GetWeight(LC::Math::WeightTag::z)->data;

    _tplane.series.GenerateDifferentials(data->directors.get(), _tplane.neighbors.get(), data->neighbors.get(), dx, dy, dz, _tplane.numNodes, data->subnodes, data->nodes, data->knn);

    _tplane.usable = true;
#else
    _plane.interpolant.ComputeFactorization(data->position.get(), _plane.neighbors.get(), *(data->RBF.get()), metric, _plane.numNodes, data->nodes, _plane.knn);
    _plane.generated = true;
#endif
}

void Sandbox::updateInterpolant() {

#if TAYLOR_METHOD
    if (!_tplane.usable) generateInterpolant();
#else
    if (!_plane.generated) generateInterpolant();
    Dataset* data = (Dataset*)_solver->GetDataPtr();
    _plane.interpolant.ComputeWeights(data->directors.get(), _plane.neighbors.get(), _plane.numNodes, data->nodes, _plane.knn);
    _plane.usable = true;
#endif
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