#include <lclab2.h>
#include <complex>
#include "Widget.h"

#define USE_PLANE 1


using namespace Magnum;
using namespace Math::Literals;

/*
    TODO:
    -   GPU compatibility
    -   GPU accelerated relax
    -   GPU accelerated POM (think computing and
        combining smaller stacks)
    -   More relax features
    -   Save file system
    -   Load file system
    -   Isosurfaces
    -   Interpolation
    -   RBF features
    -   Streamlines
    -   Q-tensor
    -   Dynamic evolution methods
    -   Steady state evaluation methods
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
    void updateColor();
    void POM();

    // Tested geometries
    LC::SphereArray _grid;
    LC::Torus _sheet;
    LC::DynamicColorSheet _dsheet;
    // Torus with PhongGL shader
    LC::NormalTorus _sheetNormal;

    /*
        GUI Widget
        *** Make sure to keep decoupled from actual simulation core functionality
    */
    Widget _widget;

    // Solver
    LC::Solver* _solver;

    LC::Header _header;

    // test save data
    char data1 = 'g';
    int data2 = 5;
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

    
    //_sheet.Init();
    //_sheetNormal.Init();


    _solver = new LC::FrankOseen::ElasticOnly::FOFDSolver;

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;
    using T4 = LC::FrankOseen::ElasticOnly::FOFDSolver::Tensor4;

    /* Setup data */
    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    data->voxels[0] = 32;
    data->voxels[1] = 32;
    data->voxels[2] = 32;

    // hard
    data->bc[0] = 0;
    data->bc[1] = 0;
    data->bc[2] = 0;

    // toron stable dimensions for one const algebraic
    data->cell_dims[0] = 0.35;
    data->cell_dims[1] = 0.35;
    data->cell_dims[2] = 0.35;

    data->k11 = LC::FrankOseen::ElasticConstants::_5CB("k11");
    data->k22 = LC::FrankOseen::ElasticConstants::_5CB("k22");
    data->k33 = LC::FrankOseen::ElasticConstants::_5CB("k33");

    // Generate a toron
    // Make a director modifying config that takes a function pointer to a lambda that modifies that director

    auto toron = [](T4 &n, int i, int j, int k, int *voxels) {
    
        int d[3] = { voxels[0] / 4, voxels[1] / 4, voxels[2] / 4 };

        if (abs(k - voxels[2] / 2) < d[2] && abs(i - voxels[0] / 2) < d[0] && abs(j - voxels[1] / 2) < d[1]) {
            n(i, j, k, 2) = -1.0;
        }
        else {
            n(i, j, k, 2) = 1.0;
        }

        n(i, j, k, 0) = 0.0;
        n(i, j, k, 1) = 0.0;
        
    };

    data->config = toron;

    _solver->Init();

    /* Setup spheres */
    _grid.NX = data->voxels[0];
    _grid.NY = data->voxels[1];
    _grid.CX = data->cell_dims[0];
    _grid.CY = data->cell_dims[1];

    _grid.Init();

    _dsheet.NX = data->voxels[0];
    _dsheet.NY = data->voxels[1];
    _dsheet.CX = data->cell_dims[0];
    _dsheet.CY = data->cell_dims[1];
    _dsheet.Init();


    // Colors
    updateColor();

    // Write data header
    // HeaderObject initializer list is { 'variable', 'size', 'location' }
    LC::Header::HeaderObject obj[] = { { "char", sizeof(char), 0 }, { "int", sizeof(int), 1 } };

    _header.headerObjects.push_back({ obj[0], 0 });
    _header.headerObjects.push_back({ obj[1], 0 });




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

    {
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::Begin("lclab2", 0, window_flags);

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("File")) {

                if (ImGui::MenuItem("New", "Ctrl+N")) {}
                if (ImGui::MenuItem("Open", "Ctrl+O")) {
                
                    auto of = pfd::open_file("Select save file name", LCLAB2_ROOT_PATH,
                        { "Text Files (.txt .text)", "*.txt *.text",
                          "All Files", "*" },
                        pfd::opt::none);

                    auto res = of.result();

                    if (!res.empty()) {
                        _header.read(res[0]);
                        // See if successful

                        for (const auto &obj : _header.headerObjects)
                            std::cout << obj.first.variable << ","
                                << obj.first.size_in_bytes << ","
                                << obj.first.location << std::endl;
                    }


                    // Read some test data
                    _header.readBody();

                    char* d1 = 0;
                    int* d2 = 0;

                    // After reading data, always keep a handle on it.
                    // Header only makes the data when it reads.
                    // It is up to the user to handle after that
                    d1 = (char*)_header.passData(0);
                    d2 = (int*)_header.passData(1);

                    std::cout << *d1 << std::endl;
                    std::cout << *d2 << std::endl;

                    delete d1;
                    delete d2;

                }

                if (ImGui::MenuItem("Save", "Ctrl+S")) {

                    // File save
                    auto sf = pfd::save_file("Select save file name", LCLAB2_ROOT_PATH,
                        { "Text Files (.txt .text)", "*.txt *.text",
                          "All Files", "*" },
                        pfd::opt::none);


                    // Pass file to header to be written
                    auto res = sf.result();

                    if (!res.empty()) {

                        // Bind data
                        _header.headerObjects[0].second = &data1;
                        _header.headerObjects[1].second = &data2;


                        _header.write(res);
                        _header.writeBody();
                    }


                }
                if (ImGui::MenuItem("Save As..")) {}

                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        if (ImGui::Button("Test Window"))
            _widget.showDemoWindow ^= true;
        if (ImGui::Button("Another Window"))
            _widget.showAnotherWindow ^= true;

        // Dropdown menu for lc types
        {
            const char* items[] = { "5CB" };
            std::string currentItem;
            if (_widget.lcType == LC::FrankOseen::LC_TYPE::_5CB) currentItem = "5CB";
            else currentItem = "Error";

            if (ImGui::BeginCombo("LC Type", currentItem.c_str())) {

                for (int n = 0; n < IM_ARRAYSIZE(items); n++) {
                    bool selected = (currentItem == items[n]);
                    if (ImGui::Selectable(items[n], selected)) {
                        _widget.lcType = static_cast<LC::FrankOseen::LC_TYPE>(n);
                        currentItem = items[n];
                    }

                    if (selected)
                        ImGui::SetItemDefaultFocus();
                }

                ImGui::EndCombo();
            }

        }

        {
            Float pitch = _widget.pitch.first;
            ImGui::InputFloat("Pitch (um)", &pitch);
            _widget.pitch.first = pitch;




            if (ImGui::CollapsingHeader("POM Settings")) {

                ImGui::Checkbox("Enable POM", &_widget.POM);
                // Dropdown menu for waveplates
                {
                    const char* items[] = { "None", "530 nm Full" };
                    std::string currentItem;
                    if (_widget.waveplate == Widget::Waveplate::Full530nm) currentItem = "530 nm Full";
                    else currentItem = "None";

                    if (ImGui::BeginCombo("Waveplate", currentItem.c_str())) {

                        for (int n = 0; n < IM_ARRAYSIZE(items); n++) {
                            bool selected = (currentItem == items[n]);
                            if (ImGui::Selectable(items[n], selected)) {
                                _widget.waveplate = static_cast<Widget::Waveplate>(n);
                                currentItem = items[n];
                            }

                            if (selected)
                                ImGui::SetItemDefaultFocus();
                        }

                        ImGui::EndCombo();
                    }

                }


                // Show additional options in a dropdown (TODO)
                std::array<Float, 3> lampIntensity = _widget.intensity;
                ImGui::InputFloat3("Lamp intensity", &lampIntensity[0]);
                _widget.intensity = lampIntensity;

                std::array<Float, 3> rgbColors = _widget.rgbColors;
                ImGui::InputFloat3("RGB", &rgbColors[0]);
                _widget.rgbColors = rgbColors;

                ImGui::Checkbox("Crossed polarizer", &_widget.crossedPolarizer);
                
                if (_widget.crossedPolarizer) {
                    Float polarizerAngle = _widget.polarizerAngle;
                    ImGui::InputFloat("Polarizer angle", &polarizerAngle);
                    _widget.polarizerAngle = polarizerAngle;
                }
            
            }
        }

        _widget.updateImage = ImGui::Button("Update Image");


        // Set Cycle
        ImGui::InputInt("Cycle", &_widget.cycle);

        // Relaxation rate
        {
            using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;

            Dataset* data = (Dataset*)(_solver->GetDataPtr());
            Float relaxRate = data->rate;
            ImGui::InputFloat("Relax rate", &relaxRate);
            data->rate = relaxRate;
        }

        // Pressed the relax button
        _widget.relax = ImGui::Button("Relax");
        ImGui::SameLine();
        _widget.print = ImGui::Button("Print");
        
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
            1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
    
        ImGui::End();
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

#if USE_PLANE
    _dsheet.Draw(_arcballCamera, _projectionMatrix);
#else
    _grid.Draw(_arcballCamera, _projectionMatrix);
#endif
    //_sheet.Draw(_arcballCamera, _projectionMatrix);
    //_sheetNormal.Draw(_arcballCamera, _projectionMatrix);

    // Make sure to draw gui last, otherwise the graphics will write over the GUI
    {
        /* Set appropriate states. If you only draw ImGui, it is sufficient to
           just enable blending and scissor test in the constructor. */
        guiRenderer();

        _imgui.drawFrame();

    }

    swapBuffers();

    if (_widget.relax) {

        _solver->Relax(_widget.cycle);
        _widget.updateImage = true;
    }
    
    if (_widget.updateImage) {
        // For now just xy midplane
        if(_widget.POM) POM();
        else updateColor();
    }
    
    if (_widget.print) {
        _solver->Print();
    }


    if (moving || _ioUpdate) redraw();
}

// Todo: inherit all these events through Application
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
    if (_imgui.handleTextInputEvent(event)) _ioUpdate = true;
}

void Sandbox::mouseReleaseEvent(MouseEvent& event) {

    if (_imgui.handleMouseReleaseEvent(event)) _ioUpdate = true;
}

void Sandbox::keyPressEvent(KeyEvent& event) {
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::keyReleaseEvent(KeyEvent& event) {
    if (_imgui.handleKeyReleaseEvent(event)) _ioUpdate = true;
}

void Sandbox::updateColor() {

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;

    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    // Colors
    std::size_t slice = data->voxels[0];
    std::size_t cross_slice = data->voxels[1] * slice;
    std::size_t volslice = data->voxels[2] * cross_slice;

    // Convenient indexing lambdas
    auto global_matlab_idx = [volslice, cross_slice, slice](int i, int j, int k, int l) {

        return volslice * l + cross_slice * k + slice * j + i;
    };

    auto cross_idx = [slice](int i, int j) {
        return slice * j + i;
    };

    LC::scalar theta, phi, nx, ny, nz;
    for (int i = 0; i < data->voxels[0]; i++) {
        for (int j = 0; j < data->voxels[1]; j++) {

            nx = data->directors[global_matlab_idx(i, j, data->voxels[2] / 2, 0)];
            ny = data->directors[global_matlab_idx(i, j, data->voxels[2] / 2, 1)];
            nz = data->directors[global_matlab_idx(i, j, data->voxels[2] / 2, 2)];
            // Compute theta and phi

            theta = acos(nz);
            phi = M_PI / 2.0f - atan2(ny, nx);

            // Compute color

            Color3 hsv;
            hsv[0] = phi / (2.0f * M_PI);
            hsv[1] = theta / M_PI * 2.0f;
            hsv[2] = 2.0f - theta / M_PI * 2.0f;

            if (hsv[0] > 1.0f) hsv[0] = 1.0f;
            if (hsv[1] > 1.0f) hsv[1] = 1.0f;
            if (hsv[2] > 1.0f) hsv[2] = 1.0f;

            _grid.sphereInstanceData[cross_idx(i, j)].color = Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] });
            _dsheet.data[cross_idx(i, j)].color = Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] });
        }
    }

    // Update sheet
    _dsheet.sheetBuffer.setData(_dsheet.data, GL::BufferUsage::DynamicDraw);
}

// Only for 5CB.
void Sandbox::POM() {

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;

    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    // Colors
    std::size_t slice = data->voxels[0];
    std::size_t cross_slice = data->voxels[1] * slice;
    std::size_t volslice = data->voxels[2] * cross_slice;

    LC::SIscalar pitch = _widget.pitch;
    LC::scalar dop = data->cell_dims[2];

    Widget::Waveplate waveplate = _widget.waveplate;

    LC::scalar thickness = pitch.first * dop;
    LC::scalar dz = thickness * 1e-6 / data->voxels[2]; // meters
    LC::scalar lambda[3]; // rgb
    LC::scalar intensity[3]; // lamp intensity

    for (int i = 0; i < 3; i++) {
        lambda[i] = _widget.rgbColors[i] * 1e-9;
        intensity[i] = _widget.intensity[i];
    }

    LC::scalar gamma = _widget.gamma;

    // Jones matrix
    Eigen::Matrix2cd M, m;

    // 5CB
    LC::scalar n0;
    LC::scalar ne;

    if (_widget.lcType == LC::FrankOseen::LC_TYPE::_5CB) {
        n0 = LC::FrankOseen::OpticalConstants::_5CB("n_o").first;
        ne = LC::FrankOseen::OpticalConstants::_5CB("n_e").first;
    }


    // Crossed polarizer (1), uncrossed (0)
    LC::scalar crossed = (LC::scalar)_widget.crossedPolarizer;
    LC::scalar polarizerAngle = _widget.polarizerAngle; // deg

    LC::scalar th = crossed * polarizerAngle * M_PI / 180.0;

    const std::complex<LC::scalar> ii(0, 1);
    
    // Convenient indexing lambdas
    auto global_matlab_idx = [volslice, cross_slice, slice](int i, int j, int k, int l) {

        return volslice * l + cross_slice * k + slice * j + i;
    };

    auto cross_idx = [slice](int i, int j) {
        return slice * j + i;
    };

    LC::scalar theta, phi, nx, ny, nz;
    for (int i = 0; i < data->voxels[0]; i++) {
        for (int j = 0; j < data->voxels[1]; j++) {

            // Polarization states for rgb colors
            Eigen::Vector2cd Eo[] = { {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0} };

            for (int k = 0; k < data->voxels[2]; k++) {

                nx = data->directors[global_matlab_idx(i, j, k, 0)];
                ny = data->directors[global_matlab_idx(i, j, k, 1)];
                nz = data->directors[global_matlab_idx(i, j, k, 2)];

                // Compute theta and phi
                theta = M_PI / 2.0 - atan2(nz, sqrt(nx * nx + ny * ny));
                phi = atan2(ny, nx);

                const LC::scalar ct = cos(theta);
                const LC::scalar st = sin(theta);
                const LC::scalar cp = cos(phi);
                const LC::scalar sp = sin(phi);


                for (int rgb = 0; rgb < 3; rgb++) {
                    const LC::scalar delta0 = 2.0 * M_PI / lambda[rgb] * dz * n0;
                    const LC::scalar ne_th = ne * n0 / sqrt(pow(n0 * st, 2.0) + pow(ne * ct, 2.0));
                    const LC::scalar deltaE = 2.0 * M_PI / lambda[rgb] * dz * ne_th;

                    const std::complex<LC::scalar> eidE = std::exp(ii * deltaE);
                    const std::complex<LC::scalar> eid0 = std::exp(ii * delta0);

                    M(0, 0) = cp * cp * eidE + sp * sp * eid0;
                    M(0, 1) = sp * cp * (eidE - eid0);
                    M(1, 0) = M(0, 1);
                    M(1, 1) = sp * sp * eidE + cp * cp * eid0;

                    Eo[rgb] = M * Eo[rgb];
                }


            }

            // Additional waveplate
            if (waveplate == Widget::Waveplate::Full530nm) {

                const std::array<int, 2> rb = { 0, 2 };
                const LC::scalar thick = 5.8889e-05 * 7.0;


                for (const auto& col : rb) {
                    const LC::scalar de_wp = 2.0 * M_PI / lambda[col] * thick * 1.55338;
                    const LC::scalar do_ep = 2.0 * M_PI / lambda[col] * thick * 1.54425;

                    m(0, 0) = 0.5 * (exp(ii * de_wp) + exp(ii * do_ep));
                    m(0, 1) = 0.5 * (exp(ii * de_wp) - exp(ii * do_ep));
                    m(1, 0) = m(0, 1);
                    m(1, 1) = m(0, 0);

                    Eo[col] = m * Eo[col];
                }

            }

            Color3 color;

            for (int rgb = 0; rgb < 3; rgb++) {
                LC::scalar I = std::abs(Eo[rgb](0) * cos(th) + Eo[rgb](1) * sin(th));
                I *= I;
                color[rgb] = pow(intensity[rgb] * I, gamma);
            }

            _grid.sphereInstanceData[cross_idx(i, j)].color = color;
            _dsheet.data[cross_idx(i, j)].color = color;
        }
    }

    // Update sheet
    _dsheet.sheetBuffer.setData(_dsheet.data, GL::BufferUsage::DynamicDraw);
}

LC::Application* LC::createApplication(int argc, char **argv) {

	return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}