#include <lclab2.h>
#include "Widget.h"

using namespace Magnum;
using namespace Math::Literals;
using FOFDSolver = LC::FrankOseen::Electric::FOFDSolver;
using Dataset = FOFDSolver::Dataset;

#define MAX_GRAPH_POINTS 1000
// Current use case for interpolating preimages is inefficient,
// use to interpolate director field THEN copy before passing to isoGenerator
#define TESTINTERP 0

enum class Axis { x = 0, y = 1, z = 2 };
struct CrossX {
    Float axisPosition;
    Axis axis;
    std::pair<Containers::Optional<GL::Mesh>, std::unique_ptr<LC::DynamicColorSheet>> section;
    bool draw = true;
    static LC::DynamicColorSheet::PositionFunction Position[3];
};
// Permutations of color sheets
LC::DynamicColorSheet::PositionFunction CrossX::Position[3] = { [](Float x, Float y, Float z) { return Vector3{ z, x, y }; },
        [](Float x, Float y, Float z) { return Vector3{ x, z, y }; },
        [](Float x, Float y, Float z) { return Vector3{ x, y, z }; } };

struct BoxInstanceData {
    Matrix4 transformationMatrix;
    Color3 color;
};

struct Preimage {
    Preimage() {}
    Preimage(float t, float p, float iso = 0.07f) : theta(t), phi(p), isoLevel(iso) {}
    
    friend bool operator == (const Preimage& p1, const Preimage& p2) {
        return (p1.theta == p2.theta && p1.phi == p2.phi) ? 1 : 0;
    }

    float theta = 0.0f;
    float phi = 0.0f;
    bool draw = true;
    float isoLevel = 0.07f;
    float alpha = 1.0f;
    bool opentab = true;
    LC::Surface surface;

    Containers::Optional<GL::Mesh> mesh;

};

class Sandbox : public LC::Application
{
public:

    explicit Sandbox(const Arguments& arguments);
    ~Sandbox();

private:
    virtual void drawEvent() override;
    void keyPressEvent(KeyEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;
    void updateColor();
    void POM();
    void initVisuals();
    void updateBoxes(bool first = false);
    void drawBoxes();
    void imageMenu();
    void handlePOMWindow();
    void handlePreimageWindow();
    void handleNonlinearImagingWindow();
    void handleLCINFOWindow();
    void handleModificationWindow();
    void generateIsosurface();
    void computeEnergy();

    template <typename T>
    void dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map);

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader;
    Shaders::PhongGL _phongShader;

    SceneGraph::DrawableGroup3D _transparentDrawables;
    SceneGraph::DrawableGroup3D _transparentNormalDrawables;

    std::unique_ptr<LC::Drawable::Object3D> _preimageManipulator;

    Containers::Array<CrossX> _crossSections;

    GL::Mesh _boxMesh{ NoCreate };
    GL::Buffer _boxInstanceBuffer{ NoCreate };
    Shaders::FlatGL3D _boxShader{ NoCreate };
    Containers::Array<BoxInstanceData> _boxInstanceData;

    /*
        GUI Widget
        *** Make sure to keep decoupled from actual simulation core functionality
    */
    Widget _widget;
    LC::Header _header;
    LC::Imaging::UniformGrid::POM _pomImager;
#if TESTINTERP
    LC::Math::Isosurface<LC::Math::Interp3Map<float>, float> _isoGenerator;
#else
    LC::Math::Isosurface<float*, float> _isoGenerator;
#endif
    std::list<Preimage> _preimages;

    LC::Imaging::ImageSeries _image_series, _nonlin_image_series;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("FO Electric Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable)

} {

    Utility::Arguments args;
    args.addOption('q', "topological-charge", "1")
        .setHelp("topological-charge", "topological charge", "Q")
        .addOption('n', "nodes-per-pitch", "30")
        .setHelp("nodes-per-pitch", "Nodes per pitch", "N")
        .addOption('v', "voltage", "1.0")
        .setHelp("voltage", "Voltage", "V")
        .addSkippedPrefix("magnum")
        .parse(arguments.argc, arguments.argv);

    _transparentShader = Shaders::VertexColorGL3D{};

    _phongShader = Shaders::PhongGL{ Shaders::PhongGL::Flag::VertexColor };


    /* Setup the GUI */
    setupGUI();

    /* Setup window and parameters */
    enableDepthTest();
    enableFaceCulling();

    /* Setup the camera */
    setupCamera(1.0f, CameraType::Group);

    /* Loop at 60 Hz max - setSwapInterval(1) sets maximum to monitor frame rate */
    setSwapInterval(0);
    setMinimalLoopPeriod(16);


    _solver = std::make_unique<FOFDSolver>();

    /* Setup data using function chaining */
    Dataset* data = (Dataset*)(_solver->GetDataPtr());



    //Dataset::Config plan = Dataset::Planar(2);
    //Dataset::Config heli = Dataset::Heliknoton(1, 0.6666, 1.135, { 0.0, 0.0, 0.0 }, false);


    ///Dataset::Config custom_cfg = [plan, heli](FOFDSolver::Tensor4& nn, int i, int j, int k, int* voxels) {
    ///    plan(nn, i, j, k, voxels);
    //    heli(nn, i, j, k, voxels);
    //};

    _widget.topological_charge = args.value<int>("topological-charge");
    _widget.npp = args.value<int>("nodes-per-pitch");
    _widget.voltage = args.value<LC::scalar>("voltage");
    int Q = _widget.topological_charge;
    _widget.celldims = { 2.f * Q + 1, 2.f * Q + 1, 2.f * Q + 1 };
    int npp = _widget.npp;
    float v0 = _widget.voltage;


    _widget.boundaries[0] = 1;
    _widget.boundaries[1] = 1;
    _widget.boundaries[2] = 0;


    /*
            Settings for nested heliknotons
    */
    (*data).Voxels((2 * Q + 1)* npp, (2 * Q + 1)* npp, (2 * Q + 1)* npp)
        .Boundaries(_widget.boundaries[0], _widget.boundaries[1], _widget.boundaries[2])
        .Cell(2 * Q + 1, 2 * Q + 1, 2 * Q + 1)
        .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        .Configuration(Dataset::Heliknoton(Q));

    /*
        Set voltage
    */

    (*data).vconfig = [v0](FOFDSolver::Tensor3& vv, int i, int j, int k, int* voxels) {
        vv(i, j, k) = v0 * (k / (voxels[2] - 1));
    };

    if (0)
    {
        // Test with Q = 2 for now
        Q = 2;

        LC::scalar dz = 0.5 * (2. * Q + 1.) / (2. * Q);

        /*
            Settings for merged structure
        */

        // Create heliknotons
        std::vector <Eigen::Matrix<LC::scalar, 3, 1>> translations;
        translations.push_back({ 0.0, 0.0, dz });
        translations.push_back({ 0.0, 0.0, -dz });

        Dataset::Config helis = Dataset::Heliknoton(1, translations, 0.50, 3);

        (*data).Voxels(5 * npp, 5 * npp, (2 * Q + 1) * npp)
            .Boundaries(1, 1, 0)
            .Cell(5, 5, 2 * Q + 1)
            .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
            .Configuration(helis);
    }


    //.Configuration(custom_cfg);

    _solver->Init();

    _relaxFuture.second = false;

    initVisuals();
    // Colors
    updateColor();

    _boxShader = Shaders::FlatGL3D{
                Shaders::FlatGL3D::Flag::VertexColor |
                Shaders::FlatGL3D::Flag::InstancedTransformation };
    _boxInstanceBuffer = GL::Buffer{};
    _boxMesh = MeshTools::compile(Primitives::cubeWireframe());
    _boxMesh.addVertexBufferInstanced(_boxInstanceBuffer, 1, 0,
        Shaders::FlatGL3D::TransformationMatrix{},
        Shaders::FlatGL3D::Color3{});

    updateBoxes(true);

    _widget.series_x_axis_vec.resize(MAX_GRAPH_POINTS);
    _widget.energy_series_vec.resize(MAX_GRAPH_POINTS);

    _image_series = LC::Imaging::ImageSeries((std::int32_t)data->voxels[0], (std::int32_t)data->voxels[1], "default-POM");


    LC_INFO("Created client application!");
}

Sandbox::~Sandbox() {
    LC_INFO("Destroying client application.");
}

/*
    Main simulation loop
*/
void Sandbox::drawEvent() {
    GL::Renderer::setClearColor(_widget.clearColor);
    GL::defaultFramebuffer.clear(
        GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    _imgui.newFrame();
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    FOFDSolver* solver = (FOFDSolver*)(_solver.get());
    {
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;

        if (_widget.showSettings) {

            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::Begin("lclab2", 0, window_flags);

            {
                std::function<void()> loadAction = [this]() {

                    initVisuals();
                    // Reset energy chart
                    Dataset* dataloc = (Dataset*)(_solver->GetDataPtr());
                    _widget.energy_series = std::list<LC::scalar>{};
                    _widget.series_x_axis = std::list<LC::scalar>{};
                    for (int i = 0; i < _widget.energy_series_vec.size(); i++) {
                        _widget.energy_series_vec[i] = 0.0;
                        _widget.series_x_axis_vec[i] = dataloc->numIterations;
                    }

                    // ** Only intended to be used for POM recording
                    _image_series = LC::Imaging::ImageSeries((std::int32_t)dataloc->voxels[0], (std::int32_t)dataloc->voxels[1], _header.readFile + "-POM");
                    for (int d = 0; d < 3; d++) _widget.celldims[d] = dataloc->cell_dims[d];

                };
                saveMenu(_widget.updateImageFromLoad, loadAction);
            }

            // Background color
            ImGui::Text("Background Color");
            ImGui::SameLine();
            ImGui::ColorEdit4("##RefColor", &_widget.clearColor[0], ImGuiColorEditFlags_NoInputs);

            imageMenu();

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
                _widget.showLCINFO ^= true;

            ImGui::SameLine();

            if (ImGui::Button("Modify LC"))
                _widget.showModificationWindow ^= true;

            


            // Dropdown menu for lc types
            {
                std::map<LC::FrankOseen::LC_TYPE, std::string> lcMap = LC::FrankOseen::LiquidCrystal::Map();
                dropDownMenu<LC::FrankOseen::LC_TYPE>("LC Type", data->lc_type, lcMap);
            }

            // Dropdown menu for relax method
            {
                using namespace LC::FrankOseen::ElasticOnly;
                auto map = Dataset::RelaxMap();
                dropDownMenu<Dataset::RelaxKind>("Relax method", data->relaxKind, map);
            }

            float t1 = _widget.energyErrorThreshold * 1e8f;

            ImGui::PushItemWidth(100.0f);
            ImGui::InputFloat("(10^-8) Energy Threshold", &t1);

            _widget.energyErrorThreshold = t1 * 1e-8f;

            ImGui::InputFloat("Voltage", &_widget.voltage);
            ImGui::InputInt("Voltage iterations", &_widget.voltage_iterations);

            ImGui::PopItemWidth();

            if (ImGui::Button("Update voltage")) {
                solver->SetVoltage(_widget.voltage, _widget.voltage_iterations);
                LC_INFO("Updated voltage");
            }
           
            handleModificationWindow();
            handleNonlinearImagingWindow();

            {
                ImGui::SliderFloat("Plane alpha", &_widget.alpha, 0.0f, 1.0f);
                
                ImGui::Checkbox("Enable POM", &_widget.POM);


                handlePOMWindow();
                handlePreimageWindow();

                {
                    bool prior = _widget.midplane;
                    ImGui::Checkbox("Midplane", &_widget.midplane);
                    if (!_widget.midplane)
                    {
                        // Initialize iplane
                        if (prior)
                            for (int d = 0; d < 3; d++)
                                _widget.iPlane[d] = data->voxels[d] / 2;

                        ImGui::SliderInt3("Image plane", &_widget.iPlane[0], 0,
                            *std::max_element(data->voxels.begin(), data->voxels.end()));

                        // Clamp within bounds
                        for (int d = 0; d < 3; d++) {
                            _widget.iPlane[d] = (_widget.iPlane[d] > data->voxels[d] - 1) ? data->voxels[d] - 1 : _widget.iPlane[d];
                            _widget.iPlane[d] = (_widget.iPlane[d] < 0) ? 0 : _widget.iPlane[d];
                        }
                    }
                }

                if (ImGui::Button("Plane Selector"))
                    ImGui::OpenPopup("plane_selector");

                std::map<std::pair<std::string, Axis>, bool&> planes{ {{"yz", Axis::x }, _crossSections[0].draw },
                    {{"xz", Axis::y }, _crossSections[1].draw },
                    {{"xy", Axis::z }, _crossSections[2].draw } };


                if (ImGui::BeginPopup("plane_selector")) {
                    ImGui::Text("Draw planes");
                    ImGui::Separator();
                    for (auto& p : planes) {
                        ImGui::MenuItem(p.first.first.c_str(), "", &p.second);
                    }
                    ImGui::EndPopup();
                }
            }

            ImGui::Text("Go to");
            ImGui::SameLine();
            if (ImGui::Button("xy"))
                _manipulator->setTransformation(Matrix4{});
            ImGui::SameLine();
            if (ImGui::Button("xz"))
                _manipulator->setTransformation(Matrix4::rotationZ(Rad(M_PI)) * Matrix4::rotationY(Rad(M_PI)) * Matrix4::rotationX(Rad(M_PI / 2.0f)));
            ImGui::SameLine();
            if (ImGui::Button("yz"))
                _manipulator->setTransformation(Matrix4::rotationY(Rad(M_PI / 2.0f)) * Matrix4::rotationX(Rad(-M_PI / 2.0f)));

            if (ImGui::Button("Update Image") || _widget.updateImageFromLoad)
                _widget.updateImage = true;

            _widget.updateImageFromLoad = false;

            // Set Cycle
            ImGui::PushItemWidth(100.0f);
            ImGui::InputInt("Cycle", &_widget.cycle);

            // Relaxation rate
            {
                Dataset* data = (Dataset*)(_solver->GetDataPtr());
                Float relaxRate = data->rate;
                ImGui::InputFloat("Relax rate", &relaxRate);
                ImGui::PopItemWidth();
                data->rate = relaxRate;
            }

            ImGui::Checkbox("GPU", &_widget.GPU);

            ImGui::SameLine();

            ImGui::Checkbox("Auto", &_widget.autoRelax);

            ImGui::RadioButton("Total Energy", &_widget.radioEn, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Total Energy Derivative", &_widget.radioEn, 1);

            // Pressed the relax button
            _widget.relax = ImGui::Button("Relax");
            ImGui::SameLine();
            ImGui::Checkbox("Compute Charge", &_widget.sampleCharge);

            if (_widget.relax) {
                _widget.updateImage = true;
            }

            // Compute change in free en:
            LC::scalar difference = 2.0 * _widget.energyErrorThreshold;
            if (_widget.energy_series.size() > 1) {
                auto it = _widget.energy_series.end();
                --it;
                auto it2 = _widget.energy_series.end();
                --it2;
                --it2;
                

                difference = *it - *it2;
            }

            // Triggered continuous relax
            if (_widget.relax && _widget.autoRelax)
                _widget.continuousRelax = true;
            // Disabled continuous relax
            else if (!_widget.autoRelax)
                _widget.continuousRelax = false;

            

            if (_relaxFuture.second) {
                ImGui::SameLine();
                ImGui::Text("Relaxing...");
                if (_widget.energy_series.size() > 1) {
                    ImGui::SameLine();
                    ImGui::Text("[dE (per cycle) = %.5e]", (float)difference / (float)_widget.cycle);
                }
            }

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));

            ImGui::End();
        }
    }

    handleLCINFOWindow();

    /* 3. Show the ImGui demo window. Most of the sample code is in
       ImGui::ShowDemoWindow() */
    if (_widget.showDemoWindow) {
        ImGui::SetNextWindowPos(ImVec2(100, 20), ImGuiCond_FirstUseEver);
        ImGui::ShowDemoWindow();
    }


   
    /* Update application cursor */
    _imgui.updateApplicationCursor(*this);


    /* Update camera */
    bool moving = true;

    if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::ArcBall))
        moving = _arcballCamera->updateTransformation();

    /* Reset state. Only needed if you want to draw something else with
        different state after. */

    if (_widget.showModificationWindow) {
        updateBoxes();
        drawBoxes();
    }

    polyRenderer();

    sortObjects(_transparentNormalDrawables);
    if (_widget.drawSurfaces)
        _camera->draw(_transparentNormalDrawables);

    // Sort objects to draw in correct order
    sortObjects(_transparentDrawables);
    _camera->draw(_transparentDrawables);


    //_sarray->Draw(_arcballCamera, _projectionMatrix);
    {
        guiRenderer();
        _imgui.drawFrame();

    }

    swapBuffers();

    if (_widget.relax && !_widget.continuousRelax) {

        // Try to relax asynchronously...
        _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, _widget.GPU);
        _relaxFuture.second = true;
        _widget.updateImage = true;
    }

    if (_widget.updateImage || _relaxFuture.second) {

        bool ready = checkRelax();

        if (ready) {
            if (_widget.updateImage) {
                updateColor();
                computeEnergy();
                // Compute topological charge
                if (_widget.sampleCharge) {
                    _widget.computed_topological_charge = LC::Math::ComputeHopfCharge(data->directors.get(), data->voxels);
                }
                _widget.updateImage = false;
            }
            
        }

    }

    // If continuous, update at the end when relax finishes
    if (_widget.continuousRelax && checkRelax()) {

        // Try to relax asynchronously...
        _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, _widget.GPU);
        _relaxFuture.second = true;
        _widget.updateImage = true;
    }

    if (moving || _ioUpdate) redraw();

}

void Sandbox::updateBoxes(bool first) {
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    if (!first) arrayResize(_boxInstanceData, 0);

    

    arrayAppend(_boxInstanceData, InPlaceInit,
        _camera->cameraMatrix() * _manipulator->transformationMatrix() *
        Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0], (float)data->cell_dims[1], (float)data->cell_dims[2] }), 0x00ffff_rgbf);

    float dx = float(_widget.shrink_interval_end[0] - _widget.shrink_interval_begin[0]) / (data->voxels[0] - 1);
    float dy = float(_widget.shrink_interval_end[1] - _widget.shrink_interval_begin[1]) / (data->voxels[1] - 1);
    float dz = float(_widget.shrink_interval_end[2] - _widget.shrink_interval_begin[2]) / (data->voxels[2] - 1);

    float xb = 0.5f * (_widget.shrink_interval_begin[0] + _widget.shrink_interval_end[0]-2) / (data->voxels[0] - 1);
    float yb = 0.5f * (_widget.shrink_interval_begin[1] + _widget.shrink_interval_end[1]-2) / (data->voxels[1] - 1);
    float zb = 0.5f * (_widget.shrink_interval_begin[2] + _widget.shrink_interval_end[2]-2) / (data->voxels[2] - 1);

    arrayAppend(_boxInstanceData, InPlaceInit,
        _camera->cameraMatrix() * _manipulator->transformationMatrix() *
        Matrix4::translation(Vector3{ (- 0.5f + xb) * (float)data->cell_dims[0],
            (-0.5f + yb) * (float)data->cell_dims[1],
            (-0.5f + zb) * (float)data->cell_dims[2] }) *
        Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0] * dx,
            (float)data->cell_dims[1] * dy,
            (float)data->cell_dims[2] * dz }), 0x00ffff_rgbf);
    
}

void Sandbox::drawBoxes() {
    _boxInstanceBuffer.setData(_boxInstanceData, GL::BufferUsage::DynamicDraw);
    _boxMesh.setInstanceCount(_boxInstanceData.size());
    _boxShader.setTransformationProjectionMatrix(_camera->projectionMatrix())
        .draw(_boxMesh);
}


void Sandbox::keyPressEvent(KeyEvent& event) {

    // Check if Ctrl + S or Ctrl + O is pressed
    if ((event.key() == KeyEvent::Key::S) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { save(); }
    else if ((event.key() == KeyEvent::Key::O) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) {
        _widget.updateImageFromLoad = open();
        if (_widget.updateImageFromLoad) initVisuals();
    }
    else {
        _widget.updateImageFromLoad = false;
    }

    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::keyReleaseEvent(KeyEvent& event) {
    if (_imgui.handleKeyReleaseEvent(event)) _ioUpdate = true;
}

void Sandbox::updateColor() {

    using T4 = LC::FrankOseen::ElasticOnly::FOFDSolver::Tensor4;
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    // Colors
    std::size_t slice = data->voxels[0];
    std::size_t cross_slice = data->voxels[1] * slice;
    std::size_t volslice = data->voxels[2] * cross_slice;
    T4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
    Float alpha = _widget.alpha;

    for (int id = 0; id < 3; id++) {
        Axis ax = static_cast<Axis>(_crossSections[id].axis);

        // Call POM and skip
        if (_widget.POM && ax == Axis::z) {
            POM();
            continue;
        }

        int xx = (ax == Axis::x) ? 1 : 0;
        int yy = (ax == Axis::y) ? 2 : xx + 1;

        std::size_t permutedSlice = data->voxels[xx];

        auto cross_idx = [permutedSlice](int i, int j) {
            return permutedSlice * j + i;
        };

        int hvox;
        if (_widget.midplane) hvox = data->voxels[id] / 2;
        else hvox = _widget.iPlane[id];

        LC::scalar theta, phi, nx, ny, nz;
        for (int i = 0; i < data->voxels[xx]; i++) {
            for (int j = 0; j < data->voxels[yy]; j++) {

                // Take an average
                if (ax == Axis::z) {
                    if (_widget.midplane) {
                        nx = (nn(i, j, hvox + 1, 0) + nn(i, j, hvox - 1, 0)) / 2.0;
                        ny = (nn(i, j, hvox + 1, 1) + nn(i, j, hvox - 1, 1)) / 2.0;
                        nz = (nn(i, j, hvox + 1, 2) + nn(i, j, hvox - 1, 2)) / 2.0;
                    }
                    else {
                        nx = nn(i, j, hvox, 0);
                        ny = nn(i, j, hvox, 1);
                        nz = nn(i, j, hvox, 2);
                    }
                }
                else if (ax == Axis::y) {
                    if (_widget.midplane) {
                        nx = (nn(i, hvox + 1, j, 0) + nn(i, hvox - 1, j, 0)) / 2.0;
                        ny = (nn(i, hvox + 1, j, 1) + nn(i, hvox - 1, j, 1)) / 2.0;
                        nz = (nn(i, hvox + 1, j, 2) + nn(i, hvox - 1, j, 2)) / 2.0;
                    }
                    else {
                        nx = nn(i, hvox, j, 0);
                        ny = nn(i, hvox, j, 1);
                        nz = nn(i, hvox, j, 2);
                    }
                }
                else if (ax == Axis::x) {
                    if (_widget.midplane) {
                        nx = (nn(hvox + 1, i, j, 0) + nn(hvox - 1, i, j, 0)) / 2.0;
                        ny = (nn(hvox + 1, i, j, 1) + nn(hvox - 1, i, j, 1)) / 2.0;
                        nz = (nn(hvox + 1, i, j, 2) + nn(hvox - 1, i, j, 2)) / 2.0;
                    }
                    else {
                        nx = nn(hvox, i, j, 0);
                        ny = nn(hvox, i, j, 1);
                        nz = nn(hvox, i, j, 2);
                    }
                }

                // NAN check

                if (nx != nx || ny != ny || nz != nz) {
                    
                    // Should check an error flag!

                    nx = 0.0;
                    ny = 0.0;
                    nz = 1.0;
                }

                // Compute theta and phi
                theta = acos(nz);
                phi = M_PI + atan2(ny, nx);

                std::size_t cidx = cross_idx(i, j);

                if (!_widget.nonlinear) _crossSections[id].section.second->vertices[cidx].color = { LC::Imaging::Colors::RungeSphere(theta, phi), alpha };
                // Green 3 photon nonlinear imaging
                else {
                    if (_widget.nonlinCircular)
                        _crossSections[id].section.second->vertices[cidx].color = { Color3::fromHsv({ Deg(120.0f), 1.0f, 1.0f - powf(nz, 6.0f) }), alpha };
                    else
                        _crossSections[id].section.second->vertices[cidx].color = { Color3::fromHsv({ Deg(120.0f), 1.0f, powf(nx * cos(M_PI / 180. * _widget.nonlinTheta) + ny * sin(M_PI / 180. * _widget.nonlinTheta), 6.0f) }), alpha };
                }
                _crossSections[id].section.second->vertices[cidx].position[id] = data->cell_dims[id] * ((float)hvox / float(data->voxels[id] - 1) - 0.5f);

                //_sarray->sphereInstanceData[cross_idx(i, j)].color = Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] });
            }
        }

        // Update sheet
        Trade::MeshData meshData = _crossSections[id].section.second->Data();
        _crossSections[id].section.second->vertexBuffer.setData(MeshTools::interleave(meshData.positions3DAsArray(), meshData.colorsAsArray()), GL::BufferUsage::DynamicDraw);



    }

}

void Sandbox::POM() {
    int i;
    for (int id = 0; id < 3; id++)
        if (_crossSections[id].axis == Axis::z)
            i = id;

    Float alpha = _widget.alpha;

    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    LC::SIscalar pitch = _widget.pitch;
    LC::scalar dop = data->cell_dims[2];

    _pomImager.n0 = LC::FrankOseen::OpticalConstants::LC(data->lc_type, LC::FrankOseen::OpticalConstants::Constant::n_o).first;
    _pomImager.ne = LC::FrankOseen::OpticalConstants::LC(data->lc_type, LC::FrankOseen::OpticalConstants::Constant::n_e).first;

    _pomImager.dop = dop;
    _pomImager.thickness = pitch.first * (dop + _pomImager.additional_layers);
    _pomImager.dz = pitch.first * dop * 1e-6 / (data->voxels[2] - 1); // meters

    _pomImager.Compute(data->directors.get(), data->voxels, (void*)(&_crossSections[i].section.second->vertices), [](void* data, const std::array<float, 4>& color, std::size_t idx) {
        Magnum::Containers::Array<LC::DynamicColorSheet::Vertex>* dataPtr = (Magnum::Containers::Array<LC::DynamicColorSheet::Vertex>*)data;
        for (int id = 0; id < 4; id++)
            (*dataPtr)[idx].color[id] = color[id];
        }, alpha);

    if (_widget.savePOM) {
        using namespace LC::Imaging;
        std::size_t sz = data->voxels[0] * data->voxels[1];
        std::unique_ptr<ImageSeries::COLOR[]> plane(new ImageSeries::COLOR[sz]);

        for (auto idx = 0; idx < sz; idx++) {
            plane[idx].R = _crossSections[i].section.second->vertices[idx].color[0] * 255;
            plane[idx].G = _crossSections[i].section.second->vertices[idx].color[1] * 255;
            plane[idx].B = _crossSections[i].section.second->vertices[idx].color[2] * 255;
            plane[idx].A = alpha * 255;
        }

        // Generate and save POM
        std::string name = getLoadFile();

        // Get just the file name
        while (1) {
            auto index1 = name.find("/");
            auto index2 = name.find("\\");

            if (index1 != std::string::npos) {
                name = name.substr(index1 + 1);
            }
            else if (index2 != std::string::npos) {
                name = name.substr(index2 + 1);
            }
            else break;

        }

        // Remove the .type if it exists
        {
            auto index = name.find(".");
            if (index != std::string::npos) {
                name = name.substr(0, index);
            }
        }

        if (name.empty()) name = "default";
        _image_series.write_file = _widget.savePOM_loc + "/" + name;
        LC_INFO("Saving POM as <{0}>", _image_series.write_file.c_str());
        _image_series.GenerateAndWriteFrame(plane.get(), data->voxels[0], data->voxels[1]);

        //auto window = Magnum::Platform::Sdl2Application::window();
        //auto windims = Magnum::Platform::Sdl2Application::windowSize();
        
        // Test writing window
        //SDL_Surface* sshot = SDL_CreateRGBSurface(0, windims.x(), windims.y(), 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);
        //SDL_Render

    }

    // Update sheet
    Trade::MeshData meshData = _crossSections[i].section.second->Data();
    _crossSections[i].section.second->vertexBuffer.setData(MeshTools::interleave(meshData.positions3DAsArray(), meshData.colorsAsArray()), GL::BufferUsage::DynamicDraw);
}


void Sandbox::initVisuals() {

    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    std::array<int, 3> voxels = data->voxels;
    std::array<LC::scalar, 3> cdims = data->cell_dims;
    _widget.shrink_interval_end = voxels;
    _widget.shrink_interval_begin = { 1, 1, 1 };

    // Create a new manipulator and set its parent
    // Reset preimageManipulator as well

    _preimageManipulator = std::make_unique<LC::Drawable::Object3D>();
    _manipulator = std::make_unique<LC::Drawable::Object3D>();
    _manipulator->setParent(&_scene);


    _preimageManipulator->setParent(_manipulator.get());

    // Set the distance from the plane
    _cameraObject.setTransformation(Matrix4::translation(Vector3::zAxis(2.2 * (std::max)(cdims[0], cdims[1]))));
    _manipulator->setTransformation(Matrix4{});
    // Manipulator rotation can be set here as well
    _crossSections = Containers::Array<CrossX>{ 3 };

    for (int id = 0; id < 3; id++) {

        Axis ax = static_cast<Axis>(id);

        int i = (ax == Axis::x) ? 1 : 0;
        int j = (ax == Axis::y) ? 2 : i + 1;

        // Start by only drawing the xy plane
        if (ax != Axis::z)
            _crossSections[id].draw = false;

        _crossSections[id].axis = ax;
        _crossSections[id].section.second = std::make_unique<LC::DynamicColorSheet>();
        _crossSections[id].section.second->Set(voxels[i], voxels[j], cdims[i], cdims[j]);
        _crossSections[id].section.second->Init(CrossX::Position[id], 0.0f);
        _crossSections[id].section.first = _crossSections[id].section.second->Mesh();
        new LC::Drawable::TransparentFlatDrawable{ *_manipulator, _transparentShader, *_crossSections[id].section.first, _crossSections[id].draw, Vector3{0.0f, 0.0f, 0.0f}, _transparentDrawables };
    }
}

void Sandbox::imageMenu() {
    if (ImGui::BeginMenuBar()) {
        if (ImGui::BeginMenu("Image")) {

            if (ImGui::MenuItem("POM Settings")) {
                // Toggle POM window
                _widget.showPOMSettings = true;
            }

            if (ImGui::MenuItem("Preimage Settings")) {
                _widget.showPreimageSettings = true;
            }

            if (ImGui::MenuItem("Nonlinear Settings")) {
                _widget.showNonlinearSettings = true;
            }

            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }
}

void Sandbox::handlePOMWindow() {

    if (_widget.showPOMSettings) {

        ImGui::Begin("Configure POM", &_widget.showPOMSettings);
        
        ImGui::Checkbox("Save POM", &_widget.savePOM);
        if (_widget.savePOM) {
            if (ImGui::Button("Choose POM location")) {
                auto res = pfd::select_folder("POM Location", LCLAB2_ROOT_PATH, pfd::opt::force_path);

                while (!res.ready(1000))
                    LC_INFO("Waiting for user...");

                if (res.ready() && !res.result().empty()) {
                    _widget.savePOM_loc = res.result();
                }
            }
            ImGui::Text("POM Location = %s", _widget.savePOM_loc.c_str());
        }

        {
            Float pitch = _widget.pitch.first;
            ImGui::InputFloat("Pitch (um)", &pitch);
            _widget.pitch.first = pitch;
        }
        // Dropdown menu for waveplates
        {
            using Waveplate = LC::Imaging::UniformGrid::POM::Waveplate;
            std::map<Waveplate, std::string> waveplateMap{ { Waveplate::None, "None" }, { Waveplate::Full530nm, "530 nm Full" } };
            dropDownMenu<Waveplate>("Waveplate", _pomImager.waveplate, waveplateMap);
        }


        // Show additional options in a dropdown (TODO)
        ImGui::InputFloat("Gamma", &_pomImager.gamma);
        ImGui::InputFloat("z-depth", &_pomImager.z_scan_ratio);
        ImGui::InputInt("+ layers (p/2)", &_pomImager.additional_layers);

        ImGui::Checkbox("Crossed polarizer", &_pomImager.polarizers);

        ImGui::InputFloat("Polarizer angle", &_pomImager.polarizerAngle);
        ImGui::End();
    }

}

void Sandbox::handlePreimageWindow() {
    if (_widget.showPreimageSettings) {
        Dataset* data = (Dataset*)(_solver->GetDataPtr());

        ImGui::Begin("Configure Preimages", &_widget.showPreimageSettings);

        ImGui::PushItemWidth(100.0f);

        if (ImGui::CollapsingHeader("Global Properties")) {

            ImGui::SliderFloat("Preimage Alpha", &_widget.preimage_alpha, 0.0f, 1.0f);
            ImGui::SameLine();
            ImGui::Checkbox("Global preimage alpha", &_widget.global_preimage_alpha);

            ImGui::PopItemWidth();
            int maxVox = (std::max)({ data->voxels[0], data->voxels[1], data->voxels[2] });
            ImGui::SliderInt3("Start Indices", &_widget.shrink_interval_begin[0], 0, maxVox);
            ImGui::SliderInt3("End Indices", &_widget.shrink_interval_end[0], 0, maxVox);
            ImGui::PushItemWidth(100.0f);

            if (!_widget.showModificationWindow) {
                updateBoxes();
                drawBoxes();
            }
        }

        ImGui::SliderInt("theta", &_widget.ptheta, 0, 180);
        ImGui::SameLine();
        ImGui::SliderInt("phi", &_widget.pphi, 0, 360);
        ImGui::SameLine();
        ImGui::SliderFloat("isoLevel", &_widget.isoLevel, 0.0f, 0.3f);


        // Clamp
        for (int d = 0; d < 3; d++) {
            if (_widget.shrink_interval_begin[d] < 1) _widget.shrink_interval_begin[d] = 1;
            if (_widget.shrink_interval_end[d] > data->voxels[d]) _widget.shrink_interval_end[d] = data->voxels[d];
            while (_widget.shrink_interval_begin[d] >= _widget.shrink_interval_end[d]) {
                --_widget.shrink_interval_begin[d];
                ++_widget.shrink_interval_end[d];
                if (_widget.shrink_interval_begin[d] < 1) _widget.shrink_interval_begin[d] = 1;
                if (_widget.shrink_interval_end[d] > data->voxels[d]) _widget.shrink_interval_end[d] = data->voxels[d];
            }
        }

        ImGui::PopItemWidth();

        if (ImGui::Button("Add Preimage")) {
            // Append a preimage to the list

            // Make sure the preimage doesn't already exist!
            bool found = false;
            for (auto& p : _preimages) {

                if (p.theta == _widget.ptheta && p.phi == _widget.pphi) {
                    found = true;
                    p.isoLevel = _widget.isoLevel;
                    
                }
            }
            if (!found)
                _preimages.emplace_back((float)_widget.ptheta, (float)_widget.pphi, _widget.isoLevel);
        }


        std::vector <std::array<float, 2>> remove_list;

        for (const auto &p : _preimages)
            if (!p.opentab) {
                remove_list.push_back({ p.theta, p.phi });
            }

        // Now we can remove
        for (const auto& p : remove_list)
            _preimages.remove(Preimage{ p[0], p[1] });
                
        // Check if using global alpha
        if (_widget.global_preimage_alpha)
            for (auto& p : _preimages) 
                p.alpha = _widget.preimage_alpha;


        // Show list to manage preimages
        ImGui::Text("(Theta, Phi) (isolevel)");
        char txtbuffer1[30];
        char txtbuffer2[30];
        if (ImGui::BeginTabBar("Selected preimages", ImGuiTabBarFlags_Reorderable)) {
            for (auto& p : _preimages) {

                int n1 = sprintf(txtbuffer1, "(%.0f, %.0f) (%f)", p.theta, p.phi, p.isoLevel);
                int n2 = sprintf(txtbuffer2, "(%.0f, %.0f)", p.theta, p.phi);
                // Copy to string
                std::string txt1, txt2;
                for (int i = 0; i < n1; i++) txt1 += txtbuffer1[i];
                for (int i = 0; i < n2; i++) txt2 += txtbuffer2[i];

                if (p.opentab && ImGui::BeginTabItem(txt1.c_str(), &p.opentab, ImGuiTabItemFlags_None))
                {
                    ImGui::Text("%s", txt1.c_str());
                    // Allow isoLevel to be adjusted
                    if (!_widget.global_preimage_alpha)
                        ImGui::SliderFloat(("Iso Alpha " + txt2).c_str(), &p.alpha, 0.0f, 1.0f);
                    ImGui::SliderFloat(("Iso Level " + txt2).c_str(), &p.isoLevel, 0.0f, 0.3f);

                    ImGui::EndTabItem();
                }
            }
            ImGui::EndTabBar();
        }


        if (ImGui::Button("Generate Isosurfaces")) {
            generateIsosurface();
        }

        ImGui::Checkbox("Draw Surfaces", &_widget.drawSurfaces);

        ImGui::End();

    }
}

void Sandbox::handleNonlinearImagingWindow() {
    if (_widget.showNonlinearSettings) {
        ImGui::Begin("Nonlinear Imaging", &_widget.showNonlinearSettings);
        ImGui::Checkbox("Nonlinear imaging", &_widget.nonlinear);

        if (ImGui::Button("Choose scan location")) {
            auto res = pfd::select_folder("Scan Location", LCLAB2_ROOT_PATH, pfd::opt::force_path);

            while (!res.ready(1000))
                LC_INFO("Waiting for user...");

            if (res.ready() && !res.result().empty()) {
                _widget.saveNonlin_loc = res.result();
            }
        }
        ImGui::Text("Scan Location = %s", _widget.saveNonlin_loc.c_str());

        if (ImGui::CollapsingHeader("Scan properties")) {
            ImGui::Text("Nonlinear scan plane");
            ImGui::RadioButton("xy", &_widget.radio_save_nonlin_xsection, 2);
            ImGui::SameLine();
            ImGui::RadioButton("yz", &_widget.radio_save_nonlin_xsection, 0);
            ImGui::SameLine();
            ImGui::RadioButton("xz", &_widget.radio_save_nonlin_xsection, 1);


            if (_widget.nonlinCircular) {
                ImGui::Text("[-] Nonlinear theta (deg)");
            }
            else {
                ImGui::PushItemWidth(100.0f);
                ImGui::SliderFloat("Nonlinear theta (deg)", &_widget.nonlinTheta, 0.0f, 365.0f);
                ImGui::PopItemWidth();
            }
            ImGui::SameLine();
            ImGui::Checkbox("Nonlinear circular polarization", &_widget.nonlinCircular);

        }
        if (ImGui::Button("Scan")) {

            // Scan the volume
            using namespace LC::FrankOseen;
            using namespace LC::Imaging;
            using T4 = ElasticOnly::FOFDSolver::Tensor4;

            Dataset* data = (Dataset*)(_solver->GetDataPtr());
            T4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
            Float alpha = _widget.alpha;
            int id = _widget.radio_save_nonlin_xsection;

            // Get file name
            std::string name = getLoadFile();

            // Get just the file name
            while (1) {
                auto index1 = name.find("/");
                auto index2 = name.find("\\");

                if (index1 != std::string::npos) {
                    name = name.substr(index1 + 1);
                }
                else if (index2 != std::string::npos) {
                    name = name.substr(index2 + 1);
                }
                else break;

            }

            // Remove the .type if it exists
            {
                auto index = name.find(".");
                if (index != std::string::npos) {
                    name = name.substr(0, index);
                }
            }

            std::string str_xsection;
            if (id == 0) str_xsection = "yz";
            else if (id == 1) str_xsection = "xz";
            else str_xsection = "xy";

            int xx = (static_cast<Axis>(id) == Axis::x) ? 1 : 0;
            int yy = (static_cast<Axis>(id) == Axis::y) ? 2 : xx + 1;

            // Determine the xsection dimensions
            int Ni = data->voxels[xx];
            int Nj = data->voxels[yy];

            // Create the image series
            _nonlin_image_series = LC::Imaging::ImageSeries((std::int32_t)Ni, (std::int32_t)Nj);

            std::string suffix;
            if (_widget.nonlinCircular) suffix = "-circ";
            else suffix = "-d" + std::to_string(_widget.nonlinTheta);

            if (name.empty()) name = "default";
            _nonlin_image_series.write_file = _widget.saveNonlin_loc + "/" + name + "-" + str_xsection + suffix;
            

            std::size_t permutedSlice = data->voxels[xx];

            auto cross_idx = [permutedSlice](int i, int j) {
                return permutedSlice * j + i;
            };

            // Determine director components
            auto dir = [id, xx, yy, nn](int i, int j, int d, int fixed) {
                float result;
                // yz plane
                if (id == 0) {
                    result = nn(fixed, i, j, d);
                }
                // xz plane
                else if (id == 1) {
                    result = nn(i, fixed, j, d);
                }
                // xy plane
                else {
                    result = nn(i, j, fixed, d);
                }
                return result;
            };



            unsigned int sz = Ni * Nj;

            std::unique_ptr<ImageSeries::COLOR[]> plane(new ImageSeries::COLOR[sz]);


            // Green 3 photon nonlinear imaging
            auto linear = [](float nx, float ny, float theta) {
                return Color3::fromHsv({ Deg(120.0f), 1.0f, powf(nx * cos(theta) + ny * sin(theta), 6.0f) });
            };
            auto circular = [](float nz) {
                return Color3::fromHsv({ Deg(120.0f), 1.0f, 1.0f - powf(nz, 6.0f) });
            };

            Color3 color;

            // Number of scans = voxels[id]
            for (int k = 0; k < data->voxels[id]; k++) {

                // Fill the plane data for each scan
                for (int i = 0; i < data->voxels[xx]; i++) {
                    for (int j = 0; j < data->voxels[yy]; j++) {
                        unsigned int idx = cross_idx(i, j);

                        // Get director components for current slice
                        float nx = dir(i, j, 0, k);
                        float ny = dir(i, j, 1, k);
                        float nz = dir(i, j, 2, k);

                        if (_widget.nonlinCircular) color = circular(nz);
                        else color = linear(nx, ny, M_PI / 180. * _widget.nonlinTheta);

                        plane[idx].R = color[0] * 255;
                        plane[idx].G = color[1] * 255;
                        plane[idx].B = color[2] * 255;
                        plane[idx].A = alpha * 255;
                    }
                }

                // Generate image

                LC_INFO("Saving nonlinear scan as <{0}>", _nonlin_image_series.write_file.c_str());
                _nonlin_image_series.GenerateAndWriteFrame(plane.get(), data->voxels[xx], data->voxels[yy]);
            }

        }
        ImGui::End();
    }
}

void Sandbox::handleLCINFOWindow() {
    if (_widget.showLCINFO) {
        using namespace LC::FrankOseen;

        Dataset* data = (Dataset*)(_solver->GetDataPtr());

        ImGui::SetNextWindowSize(ImVec2(500, 100), ImGuiCond_FirstUseEver);
        ImGui::Begin("LC Info", &_widget.showLCINFO);


        // Relaxation time series data
        {
            int points = _widget.series_x_axis.size();
            
            if (ImPlot::BeginPlot("Free Energy Density v. Iterations", "Iterations", "Free energy density")) {
                ImPlot::SetNextMarkerStyle(ImPlotMarker_None);

                // Reset points
                if (_widget.updateImageFromLoad) points = MAX_GRAPH_POINTS;
                ImPlot::PlotLine("F", &_widget.series_x_axis_vec[0], &_widget.energy_series_vec[0], points);
                ImPlot::EndPlot();
            }
        }

        std::map<LC_TYPE, std::string> lcMap = LC::FrankOseen::LiquidCrystal::Map();
        auto relaxMethodMap = Dataset::RelaxMap();

        float dV = 0.0f;

        if (data->voltage.get()) {
            FOFDSolver::Tensor3 voltage(data->voltage.get(), data->voxels[0], data->voxels[1], data->voxels[2]);
            dV = voltage(0, 0, data->voxels[2] - 1) - voltage(0, 0, 0);
        }

        // Show LC information
        ImGui::Text("LC = %s", lcMap[data->lc_type].c_str());
        ImGui::Text("Topological Charge = %.3f", _widget.computed_topological_charge);
        ImGui::Text("Voltage = %f", dV);
        ImGui::Text("Voxels = %d %d %d", data->voxels[0], data->voxels[1], data->voxels[2]);
        ImGui::Text("Cell dimensions = %.2f %.2f %.2f", data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);
        ImGui::Text("Periodic Boundaries = %d %d %d", data->bc[0], data->bc[1], data->bc[2]);
        ImGui::Text("Iterations = %d", data->numIterations);
        ImGui::Text("Relax method = %s", relaxMethodMap[data->relaxKind].c_str());
        ImGui::End();
    }
}

void Sandbox::handleModificationWindow() {
    if (_widget.showModificationWindow) {
        FOFDSolver* solver = (FOFDSolver*)(_solver.get());
        Dataset* data = (Dataset*)(_solver->GetDataPtr());
        ImGui::SetNextWindowContentSize({ 150,100 });
        ImGui::SetNextWindowPos({ 0, 450 }, ImGuiCond_FirstUseEver);
        ImGui::Begin("LC Modification window", &_widget.showModificationWindow);

        ImGui::InputFloat3("Cell dims", &_widget.celldims[0]);
        ImGui::InputInt3("PBCs", &_widget.boundaries[0]);
        ImGui::InputInt("npp", &_widget.npp);
        ImGui::InputInt("Q", &_widget.topological_charge);

        int Q = _widget.topological_charge;
        int npp = _widget.npp;
        float v0 = _widget.voltage;

        if (ImGui::Button("Generate")) {

            (*data).Voxels(_widget.celldims[0] * npp, _widget.celldims[1] * npp, _widget.celldims[2] * npp)
                .Boundaries(_widget.boundaries[0], _widget.boundaries[1], _widget.boundaries[2])
                .Cell(_widget.celldims[0], _widget.celldims[1], _widget.celldims[2])
                .Configuration(Dataset::Heliknoton(Q));

            // Reset data
            (*data).numIterations = 0;

            _widget.energy_series = std::list<LC::scalar>{};
            _widget.series_x_axis = std::list<LC::scalar>{};
            for (int i = 0; i < _widget.energy_series_vec.size(); i++) {
                _widget.energy_series_vec[i] = 0.0;
                _widget.series_x_axis_vec[i] = 0.0;
            }

            // Initialize
            _solver->Init();
            // Update visuals
            initVisuals();
            //updateColor();
            _widget.updateImageFromLoad = true;
        }

        ImGui::SameLine();

        if (ImGui::Button("Edit existing")) {
            (*data).Boundaries(_widget.boundaries[0], _widget.boundaries[1], _widget.boundaries[2])
                .Cell(_widget.celldims[0], _widget.celldims[1], _widget.celldims[2]);
            initVisuals();
            //updateColor();
            _widget.updateImageFromLoad = true;
        }

        ImGui::InputInt("Interpolate", &_widget.interpolate);
        if (ImGui::Button("Begin Interpolate") && _widget.interpolate > 1) {

            // Create a new ptr
            unsigned int numDirectors = data->voxels[0] * data->voxels[1] * data->voxels[2];
            std::unique_ptr<LC::scalar[]> field_nn(new LC::scalar[3 * numDirectors * _widget.interpolate * _widget.interpolate * _widget.interpolate]);
            std::unique_ptr<LC::scalar[]> field_vv(new LC::scalar[numDirectors * _widget.interpolate * _widget.interpolate * _widget.interpolate]);
            // Check voltage quickly
            if (!data->voltage.get()) {
                solver->SetVoltage(1.0);
            }
            LC::Math::Interp3(data->directors.get(), field_nn.get(), data->voxels, { _widget.interpolate, _widget.interpolate , _widget.interpolate }, 3);
            LC::Math::Interp3(data->voltage.get(), field_vv.get(), data->voxels, { _widget.interpolate, _widget.interpolate , _widget.interpolate }, 1);

            // Replace data
            for (int d = 0; d < 3; d++) data->voxels[d] *= _widget.interpolate;
            numDirectors = data->voxels[0] * data->voxels[1] * data->voxels[2];
            data->directors = std::unique_ptr<LC::scalar[]>(field_nn.release());
            data->voltage = std::unique_ptr<LC::scalar[]>(field_vv.release());
            // Set energy again
            data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[numDirectors]);

            solver->Normalize();
            initVisuals();
            updateColor();
        }

        if (ImGui::Button("Helical field")) {
            auto helical = LC::Math::Planar(2 * data->cell_dims[2], data->cell_dims[2]);
            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

            std::array<float, 3> pos;

            LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
            LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
            LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);

            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; i++)
                for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; j++)
                    for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; k++) {

                        pos[0] = -data->cell_dims[0] / 2.0 + i * dx;
                        pos[1] = -data->cell_dims[1] / 2.0 + j * dy;
                        pos[2] = -data->cell_dims[2] / 2.0 + k * dz;

                        auto n = helical(pos[0], pos[1], pos[2]);

                        for (int d = 0; d < 3; d++) {
                            nn(i, j, k, d) = n[d];
                        }
                    }
        }

        // Position in [-Lx/2:Lx/2, -Ly/2:Ly/2, -Lz/2:Lz/2]                

        static std::array<float, 3> fp_arr{ 0.0, 0.0, 0.0 }, f_orien{ 0.0, 0.0,-1.0 };
        static std::array<int, 3> iperm_arr{ 0, 1, 2 };
        static float f_r1{ 0.5 }, f_r2{ 0.1 };
        static float deform = 0.0;

        ImGui::InputFloat3("Position", &fp_arr[0]);
        ImGui::InputFloat3("Orientation", &f_orien[0]);
        ImGui::InputFloat("Deformation", &deform);
        ImGui::InputInt3("sigma(xyz)", &iperm_arr[0]);

        ImGui::InputFloat("R", &f_r1);
        ImGui::InputFloat("r", &f_r2);




        if (ImGui::Button("Add torus")) {

            for (int d = 0; d < 3; d++) {

                // clamp
                if (iperm_arr[d] > 2) iperm_arr[d] = 2;
                if (iperm_arr[d] < 0) iperm_arr[d] = 0;

                // No repeats
                if (iperm_arr[d] == iperm_arr[(d + 1) % 3] || iperm_arr[d] == iperm_arr[(d + 2) % 3]) {
                    // Reset
                    iperm_arr[0] = 0;
                    iperm_arr[1] = 1;
                    iperm_arr[2] = 2;
                }
            }

            std::array<LC::scalar, 3> p_arr;
            std::array<LC::scalar, 3> orientation;
            LC::scalar r1, r2;

            for (int d = 0; d < 3; d++) {
                p_arr[d] = fp_arr[d];
                orientation[d] = f_orien[d];
            }

            r1 = f_r1;
            r2 = f_r2;

            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

            std::array<float, 3> pos;

            LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
            LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
            LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);

            auto torus = LC::Math::RubberTorus(p_arr, r1, r2, deform);

            for (int i = 0; i < data->voxels[0]; i++)
                for (int j = 0; j < data->voxels[1]; j++)
                    for (int k = 0; k < data->voxels[2]; k++) {

                        pos[iperm_arr[0]] = -data->cell_dims[0] / 2.0 + i * dx;
                        pos[iperm_arr[1]] = -data->cell_dims[1] / 2.0 + j * dy;
                        pos[iperm_arr[2]] = -data->cell_dims[2] / 2.0 + k * dz;

                        // Is inside torus
                        if (torus(pos[0], pos[1], pos[2])) {
                            nn(i, j, k, 0) = orientation[0];
                            nn(i, j, k, 1) = orientation[1];
                            nn(i, j, k, 2) = orientation[2];
                        }

                    }
        }

        int maxVox = (std::max)({ data->voxels[0], data->voxels[1], data->voxels[2] });
        ImGui::SliderInt3("Start Indices", &_widget.shrink_interval_begin[0], 0, maxVox);
        ImGui::SliderInt3("End Indices", &_widget.shrink_interval_end[0], 0, maxVox);

        // Clamp
        for (int d = 0; d < 3; d++) {
            if (_widget.shrink_interval_begin[d] < 1) _widget.shrink_interval_begin[d] = 1;
            if (_widget.shrink_interval_end[d] > data->voxels[d]) _widget.shrink_interval_end[d] = data->voxels[d];
            while (_widget.shrink_interval_begin[d] >= _widget.shrink_interval_end[d]) {
                --_widget.shrink_interval_begin[d];
                ++_widget.shrink_interval_end[d];
                if (_widget.shrink_interval_begin[d] < 1) _widget.shrink_interval_begin[d] = 1;
                if (_widget.shrink_interval_end[d] > data->voxels[d]) _widget.shrink_interval_end[d] = data->voxels[d];
            }
        }

        if (ImGui::Button("Resize")) {
            // Create a new ptr
            unsigned int numDirectors = data->voxels[0] * data->voxels[1] * data->voxels[2];
            std::array<int, 3> vNew;
            for (int d = 0; d < 3; d++) {
                vNew[d] = _widget.shrink_interval_end[d] - _widget.shrink_interval_begin[d] + 1;
                // Update cell dims
                data->cell_dims[d] *= (LC::scalar)vNew[d] / (LC::scalar)data->voxels[d];
            }




            std::unique_ptr<LC::scalar[]> field_nn(new LC::scalar[3 * vNew[0] * vNew[1] * vNew[2]]);
            std::unique_ptr<LC::scalar[]> field_vv(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);

            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
            // Check if voltage is initialized
            if (!data->voltage.get()) {
                solver->SetVoltage(1.0);
            }
            FOFDSolver::Tensor3 vv(data->voltage.get(), data->voxels[0], data->voxels[1], data->voxels[2]);

            FOFDSolver::Tensor4 nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);
            FOFDSolver::Tensor3 vv_new(field_vv.get(), vNew[0], vNew[1], vNew[2]);

            for (int i = _widget.shrink_interval_begin[0]-1; i < _widget.shrink_interval_end[0]; i++)
                for (int j = _widget.shrink_interval_begin[1]-1; j < _widget.shrink_interval_end[1]; j++)
                    for (int k = _widget.shrink_interval_begin[2]-1; k < _widget.shrink_interval_end[2]; k++) {


                        int ip = i - _widget.shrink_interval_begin[0]+1;
                        int jp = j - _widget.shrink_interval_begin[1]+1;
                        int kp = k - _widget.shrink_interval_begin[2]+1;

                        for (int d = 0; d < 3; d++)
                            nn_new(ip, jp, kp, d) = nn(i, j, k, d);

                        vv_new(ip, jp, kp) = vv(i, j, k);

                    }

            // Replace data
            for (int d = 0; d < 3; d++) data->voxels[d] = vNew[d];
            data->directors = std::unique_ptr<LC::scalar[]>(field_nn.release());
            data->voltage = std::unique_ptr<LC::scalar[]>(field_vv.release());
            data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0]* vNew[1]* vNew[2]]);

            initVisuals();
            _widget.updateImageFromLoad = true;
        }

        ImGui::SameLine();

        // Duplicate the volume along the x axis
        if (ImGui::Button("Repeat (x)")) {
            auto voxi = data->voxels;
            auto vNew = voxi;
            vNew[0] *= 2;

            std::unique_ptr<LC::scalar[]> field_nn(new LC::scalar[3 * vNew[0] * vNew[1] * vNew[2]]);
            std::unique_ptr<LC::scalar[]> field_vv(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);

            FOFDSolver::Tensor4 nn(data->directors.get(), voxi[0], voxi[1], voxi[2], 3);
            // Check if voltage is initialized
            if (!data->voltage.get()) {
                solver->SetVoltage(1.0);
            }
            FOFDSolver::Tensor3 vv(data->voltage.get(), data->voxels[0], data->voxels[1], data->voxels[2]);

            FOFDSolver::Tensor4 nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);
            FOFDSolver::Tensor3 vv_new(field_vv.get(), vNew[0], vNew[1], vNew[2]);

            for (int rep = 0; rep < 2; rep++)
                for (int i = 0; i < voxi[0]; i++)
                    for (int j = 0; j < voxi[1]; j++)
                        for (int k = 0; k < voxi[2]; k++) {

                            for (int d = 0; d < 3; d++)
                                nn_new(i + voxi[0] * rep, j, k, d) = nn(i, j, k, d);

                        }



            LC::scalar v0 = vv(0, 0, voxi[2] - 1);
            data->cell_dims[0] *= 2;

            // Replace data
            for (int d = 0; d < 3; d++) data->voxels[d] = vNew[d];
            data->directors = std::unique_ptr<LC::scalar[]>(field_nn.release());
            data->voltage = std::unique_ptr<LC::scalar[]>(field_vv.release());
            data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);
            solver->SetVoltage(v0);
            initVisuals();
            _widget.updateImageFromLoad = true;

        }
        ImGui::SameLine();

        // Duplicate the volume along the y axis
        if (ImGui::Button("Repeat (y)")) {
            auto voxi = data->voxels;
            auto vNew = voxi;
            vNew[1] *= 2;

            std::unique_ptr<LC::scalar[]> field_nn(new LC::scalar[3 * vNew[0] * vNew[1] * vNew[2]]);
            std::unique_ptr<LC::scalar[]> field_vv(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);

            FOFDSolver::Tensor4 nn(data->directors.get(), voxi[0], voxi[1], voxi[2], 3);
            // Check if voltage is initialized
            if (!data->voltage.get()) {
                solver->SetVoltage(1.0);
            }
            FOFDSolver::Tensor3 vv(data->voltage.get(), data->voxels[0], data->voxels[1], data->voxels[2]);

            FOFDSolver::Tensor4 nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);
            FOFDSolver::Tensor3 vv_new(field_vv.get(), vNew[0], vNew[1], vNew[2]);

            for (int rep = 0; rep < 2; rep++)
                for (int i = 0; i < voxi[0]; i++)
                    for (int j = 0; j < voxi[1]; j++)
                        for (int k = 0; k < voxi[2]; k++) {

                            for (int d = 0; d < 3; d++)
                                nn_new(i, j + voxi[1] * rep, k, d) = nn(i, j, k, d);

                        }



            LC::scalar v0 = vv(0, 0, voxi[2] - 1);
            data->cell_dims[1] *= 2;

            // Replace data
            for (int d = 0; d < 3; d++) data->voxels[d] = vNew[d];
            data->directors = std::unique_ptr<LC::scalar[]>(field_nn.release());
            data->voltage = std::unique_ptr<LC::scalar[]>(field_vv.release());
            data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);
            solver->SetVoltage(v0);
            initVisuals();
            _widget.updateImageFromLoad = true;

        }

        ImGui::SameLine();

        // Duplicate the volume along the z axis
        if (ImGui::Button("Repeat (z)")) {
            auto voxi = data->voxels;
            auto vNew = voxi;
            vNew[2] *= 2;

            std::unique_ptr<LC::scalar[]> field_nn(new LC::scalar[3 * vNew[0] * vNew[1] * vNew[2]]);
            std::unique_ptr<LC::scalar[]> field_vv(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);

            FOFDSolver::Tensor4 nn(data->directors.get(), voxi[0], voxi[1], voxi[2], 3);
            // Check if voltage is initialized
            if (!data->voltage.get()) {
                solver->SetVoltage(1.0);
            }
            FOFDSolver::Tensor3 vv(data->voltage.get(), data->voxels[0], data->voxels[1], data->voxels[2]);

            FOFDSolver::Tensor4 nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);
            FOFDSolver::Tensor3 vv_new(field_vv.get(), vNew[0], vNew[1], vNew[2]);

            for (int rep = 0; rep < 2; rep++)
                for (int i = 0; i < voxi[0]; i++)
                    for (int j = 0; j < voxi[1]; j++)
                        for (int k = 0; k < voxi[2]; k++) {

                            for (int d = 0; d < 3; d++)
                                nn_new(i, j, k + voxi[2] * rep, d) = nn(i, j, k, d);

                        }

            

            LC::scalar v0 = 2. * vv(0, 0, voxi[2] - 1);
            data->cell_dims[2] *= 2;

            // Replace data
            for (int d = 0; d < 3; d++) data->voxels[d] = vNew[d];
            data->directors = std::unique_ptr<LC::scalar[]>(field_nn.release());
            data->voltage = std::unique_ptr<LC::scalar[]>(field_vv.release());
            data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);
            solver->SetVoltage(v0);
            initVisuals();
            _widget.updateImageFromLoad = true;

        }



        ImGui::End();
    }

}

void Sandbox::generateIsosurface() {
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
#if TESTINTERP
    int nterp = 2;
#endif
    std::array<LC::scalar, 3> cell =
#if TESTINTERP
    { data->cell_dims[0] / (data->voxels[0] * nterp - 1),data->cell_dims[1] / (data->voxels[1] * nterp - 1),data->cell_dims[2] / (data->voxels[2] * nterp - 1) };
#else
    { data->cell_dims[0] / (data->voxels[0] - 1), data->cell_dims[1] / (data->voxels[1] - 1), data->cell_dims[2] / (data->voxels[2] - 1) };
#endif
    std::array<LC::scalar, 3> N;

    // ==========================
    unsigned int numDirectors = data->voxels[0] * data->voxels[1] * data->voxels[2];
    std::array<int, 3> vNew;
    for (int d = 0; d < 3; d++)
        vNew[d] = _widget.shrink_interval_end[d] - _widget.shrink_interval_begin[d] + 1;



    unsigned int size = vNew[0] * vNew[1] * vNew[2];
    std::unique_ptr<LC::scalar[]> field_nn(new LC::scalar[3 * size]);
    FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
    FOFDSolver::Tensor4 nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);

    for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; i++)
        for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; j++)
            for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; k++) {


                int ip = i - _widget.shrink_interval_begin[0] + 1;
                int jp = j - _widget.shrink_interval_begin[1] + 1;
                int kp = k - _widget.shrink_interval_begin[2] + 1;

                for (int d = 0; d < 3; d++)
                    nn_new(ip, jp, kp, d) = nn(i, j, k, d);

            }

    std::unique_ptr<float[]> field(new float[size]);


    for (int d = 0; d < 3; d++)
#if TESTINTERP
        _widget.preimage_translate[d] = 0.5f * (_widget.shrink_interval_end[d] - data->voxels[d] + _widget.shrink_interval_begin[d] - 1.0f) * cell[d];
#else
        _widget.preimage_translate[d] = 0.5f * (_widget.shrink_interval_end[d] - data->voxels[d] + _widget.shrink_interval_begin[d] - 1.0f) * cell[d];
#endif
    
    // ======================================================
    // Reset the manipulator
    _preimageManipulator = std::make_unique<LC::Drawable::Object3D>();
    _preimageManipulator->setParent(_manipulator.get());


    for (Preimage& pimage : _preimages) {
        // Compute the field (diffmag)
        float theta = pimage.theta / 180.f * M_PI;
        float phi = pimage.phi / 180.f * M_PI;

        N = { sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
        for (unsigned int i = 0; i < size; i++) {
            field[i] = (field_nn[i] - N[0]) * (field_nn[i] - N[0]) +
                (field_nn[i + size] - N[1]) * (field_nn[i + size] - N[1]) +
                (field_nn[i + 2 * size] - N[2]) * (field_nn[i + 2 * size] - N[2]);
        }

        std::array<float, 4> color;
        {
            auto col = LC::Imaging::Colors::RungeSphere(theta, phi);
            for (int d = 0; d < 3; d++)
                color[d] = col[d];

            color[3] = pimage.alpha;
        }

#if TESTINTERP
        // First pass field to interpolator...
        LC::Math::Interp3Map<float> fieldInt(field.get(), vNew, { nterp,nterp,nterp });
        for (int d = 0; d < 3; d++) vNew[d] *= nterp;

        _isoGenerator.GenerateSurface(fieldInt, pimage.isoLevel, vNew, cell, color);
#else
        _isoGenerator.GenerateSurface(field.get(), pimage.isoLevel, vNew, cell, color);
#endif

        if (_isoGenerator.isSurfaceValid()) {

            unsigned int nVert = _isoGenerator.NumSurfaceVertices();
            unsigned int nInd = _isoGenerator.NumSurfaceIndices();

            LC::Math::IsoVertex* verts = _isoGenerator.ReleaseSurfaceVertices();
            unsigned int* indices = _isoGenerator.ReleaseSurfaceIndices();

            LC_INFO("Successfully generated surface (verts = {0}, indices = {1})", nVert, nInd);

            // Fill magnum class with generated surface data
            pimage.surface.Init((LC::Surface::Vertex*)verts, nVert, indices, nInd, _widget.preimage_translate);

            // Delete data
            delete[] verts;
            delete[] indices;

            // Reset generator
            _isoGenerator.DeleteSurface();

            pimage.mesh = pimage.surface.Mesh();

            // Add mesh to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_preimageManipulator, _phongShader, *pimage.mesh, pimage.draw, _transparentNormalDrawables };
        }

    }
}

void Sandbox::computeEnergy() {
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    FOFDSolver* solver = (FOFDSolver*)(_solver.get());

    // Update list
    LC::scalar unit_density = data->cell_dims[0] * data->cell_dims[1] * data->cell_dims[2] / ((data->voxels[0] - 1) * (data->voxels[1] - 1) * (data->voxels[2] - 1));

    LC::scalar energy = 0.0;

    if (_widget.radioEn == 0)
        energy = solver->TotalEnergy();
    else if (_widget.radioEn == 1)
        energy = solver->TotalEnergyFunctionalDerivativeAbsSum();

    _widget.energy_series.push_back(energy * unit_density);
    _widget.series_x_axis.push_back(data->numIterations);

    // Update time series for free energy

    if (_widget.energy_series.size() == MAX_GRAPH_POINTS) {
        _widget.energy_series.pop_front();
        _widget.series_x_axis.pop_front();
    }

    // Feed lists into vectors
    int i = 0;
    for (const auto& en : _widget.energy_series) {
        _widget.energy_series_vec[i++] = en;
    }

    // Repeat for x-axis
    i = 0;
    for (const auto& it : _widget.series_x_axis) {
        _widget.series_x_axis_vec[i++] = it;
    }
}

// Menu template to make menu creation a breeze with maps
// Should really be in a separate header file in Sandbox to clear up
// the main application file...
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