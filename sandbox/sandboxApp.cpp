#include <lclab2.h>
#include "Widget.h"
using namespace Magnum;
using namespace Math::Literals;
using FOFDSolver = LC::FrankOseen::ElasticOnly::FOFDSolver;
using Dataset = FOFDSolver::Dataset;

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

    template <typename T>
    void dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map);

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader;
    SceneGraph::DrawableGroup3D _transparentDrawables;

    Containers::Array<CrossX> _crossSections;

    /*
        GUI Widget
        *** Make sure to keep decoupled from actual simulation core functionality
    */
    Widget _widget;

    LC::Header _header;
    LC::Imaging::UniformGrid::POM _pomImager;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("FOFD Elastic Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable)

} {
   
    Utility::Arguments args;
    args.addOption('q', "topological-charge", "1")
        .setHelp("topological-charge", "topological charge", "Q")
        .addOption('n', "nodes-per-pitch", "30")
        .setHelp("nodes-per-pitch", "Nodes per pitch", "N")
        .addSkippedPrefix("magnum")
        .parse(arguments.argc, arguments.argv);

    _transparentShader = Shaders::VertexColorGL3D{};

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

    int Q = args.value<int>("topological-charge");
    int npp = args.value<int>("nodes-per-pitch");


    /*
            Settings for nested heliknotons
    */
    (*data).Voxels((2 * Q + 1) * npp, (2 * Q + 1) * npp, (2 * Q + 1) * npp)
        .Boundaries(1, 1, 0)
        .Cell(2 * Q + 1, 2 * Q + 1, 2 * Q + 1)
        .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        .Configuration(Dataset::Heliknoton(Q));

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

    LC_INFO("Created client application!");
}

Sandbox::~Sandbox() {
    LC_INFO("Destroying client application.");
}

/*
    Main simulation loop
*/
void Sandbox::drawEvent() {
    GL::defaultFramebuffer.clear(
        GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    _imgui.newFrame();

    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    {
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;


        if (_widget.showSettings) {

            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::Begin("lclab2", 0, window_flags);

            {
                std::function<void()> p = [this]() {
                    initVisuals();
                };
                saveMenu(_widget.updateImageFromLoad, p);
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

            static bool p_modification_window = false;
            if (ImGui::Button("Modify LC"))
                p_modification_window ^= true;

            if (p_modification_window) {
                ImGui::SetNextWindowContentSize({ 150,100 });
                ImGui::SetNextWindowPos({ 0, 450 });
                ImGui::Begin("LC Modification window", &p_modification_window);

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

                ImGui::End();
            }
            

            {
                ImGui::SliderFloat("Plane alpha", &_widget.alpha, 0.0f, 1.0f);
                ImGui::Checkbox("Nonlinear imaging", &_widget.nonlinear);

                if (_widget.nonlinear) {
                    ImGui::SliderFloat("Nonlinear theta (deg)", &_widget.nonlinTheta, 0.0f, 365.0f);
                    ImGui::Checkbox("Nonlinear circular polarization", &_widget.nonlinCircular);
                }

                if (ImGui::CollapsingHeader("POM Settings")) {

                    ImGui::Checkbox("Enable POM", &_widget.POM);

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
                    ImGui::InputFloat3("Lamp intensity", &_pomImager.intensity[0]);
                    ImGui::InputFloat3("RGB", &_pomImager.lightRGB[0]);
                    ImGui::InputFloat("Gamma", &_pomImager.gamma);
                    ImGui::InputFloat("z-depth", &_pomImager.z_scan_ratio);
                    ImGui::InputInt("+ layers (p/2)", &_pomImager.additional_layers);

                    ImGui::Checkbox("Crossed polarizer", &_pomImager.polarizers);

                    ImGui::InputFloat("Polarizer angle", &_pomImager.polarizerAngle);
                }

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
                    for (const auto& p : planes) {
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

            _widget.updateImage = ImGui::Button("Update Image") || _widget.updateImageFromLoad;

            // Set Cycle
            ImGui::InputInt("Cycle", &_widget.cycle);

            // Relaxation rate
            {
                Dataset* data = (Dataset*)(_solver->GetDataPtr());
                Float relaxRate = data->rate;
                ImGui::InputFloat("Relax rate", &relaxRate);
                data->rate = relaxRate;
            }

            ImGui::Checkbox("GPU", &_widget.GPU);

            // Pressed the relax button
            _widget.relax = ImGui::Button("Relax");
            if (_relaxFuture.second) {
                ImGui::SameLine();
                ImGui::Text("Relaxing...");
            }

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));

            ImGui::End();
        }
    }

    /* 2. Show another simple window, now using an explicit Begin/End pair */
    if (_widget.showAnotherWindow) {
        using namespace LC::FrankOseen;

        ImGui::SetNextWindowSize(ImVec2(500, 100), ImGuiCond_FirstUseEver);
        ImGui::Begin("LC Info", &_widget.showAnotherWindow);

        std::map<LC_TYPE, std::string> lcMap = LC::FrankOseen::LiquidCrystal::Map();
        auto relaxMethodMap = ElasticOnly::FOFDSolver::Dataset::RelaxMap();

        // Show LC information
        ImGui::Text("LC = %s", lcMap[data->lc_type].c_str());
        ImGui::Text("Voxels = %d %d %d", data->voxels[0], data->voxels[1], data->voxels[2]);
        ImGui::Text("Cell dimensions = %.2f %.2f %.2f", data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);
        ImGui::Text("Periodic Boundaries = %d %d %d", data->bc[0], data->bc[1], data->bc[2]);
        ImGui::Text("Iterations = %d", data->numIterations);
        ImGui::Text("Relax method = %s", relaxMethodMap[data->relaxKind].c_str());
        ImGui::End();
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
    bool moving = true;

    if (_cameraType == CameraType::ArcBall)
        moving = _arcballCamera->updateTransformation();

    /* Reset state. Only needed if you want to draw something else with
        different state after. */

    polyRenderer();

    // Sort objects to draw in correct order
    sortObjects(_transparentDrawables);
    _camera->draw(_transparentDrawables);

    //_sarray->Draw(_arcballCamera, _projectionMatrix);
    {
        guiRenderer();
        _imgui.drawFrame();

    }

    swapBuffers();


    if (_widget.relax) {

        // Try to relax asynchronously...
        _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, _widget.GPU);
        _relaxFuture.second = true;
        _widget.updateImage = true;
    }

    if (_widget.updateImage || _relaxFuture.second) {

        bool ready = checkRelax();

        if (ready) {
            updateColor();
            _widget.updateImageFromLoad = false;
        }

    }

    if (moving || _ioUpdate) redraw();

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
                        _crossSections[id].section.second->vertices[cidx].color = { Color3::fromHsv({ Deg(120.0f), 1.0f, powf(nx * cos(M_PI / 180. * _widget.nonlinTheta)+ ny * sin(M_PI / 180. * _widget.nonlinTheta), 6.0f) }), alpha };
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
    _pomImager.dz = pitch.first * dop * 1e-6 / (data->voxels[2]-1); // meters

    _pomImager.Compute(data->directors.get(), data->voxels, (void*)(&_crossSections[i].section.second->vertices), [](void* data, const std::array<float, 4>& color, std::size_t idx) {
        Magnum::Containers::Array<LC::DynamicColorSheet::Vertex>* dataPtr = (Magnum::Containers::Array<LC::DynamicColorSheet::Vertex>*)data;
        for (int id = 0; id < 4; id++)
            (*dataPtr)[idx].color[id] = color[id];
        }, alpha);

    // Update sheet
    Trade::MeshData meshData = _crossSections[i].section.second->Data();
    _crossSections[i].section.second->vertexBuffer.setData(MeshTools::interleave(meshData.positions3DAsArray(), meshData.colorsAsArray()), GL::BufferUsage::DynamicDraw);
}


void Sandbox::initVisuals() {

    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    std::array<int, 3> voxels = data->voxels;
    std::array<LC::scalar, 3> cdims = data->cell_dims;


    // Create a new manipulator and set its parent
    _manipulator = std::make_unique<LC::Drawable::Object3D>();
    _manipulator->setParent(&_scene);
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