#include <lclab2.h>
#include <complex>
#include "Widget.h"

using namespace Magnum;
using namespace Math::Literals;

/*
    TODO:
    -   GPU compatibility
    -   GPU accelerated relax
    -   GPU accelerated POM (think computing and
        combining smaller stacks)
    -   More relax features
    -   Isosurfaces
    -   Interpolation
    -   RBF features
    -   Streamlines
    -   Q-tensor

    Ambitions:
    -   Dynamic evolution methods
    -   Steady state evaluation methods
*/

enum class Axis {
    x = 0,
    y = 1,
    z = 2
};
struct CrossX {
    Float axisPosition;
    std::unique_ptr<LC::DynamicColorSheet> section;
    bool draw = true;
    static LC::DynamicColorSheet::PositionFunction Position[3];

    void Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball,
        const Magnum::Matrix4& projection) { section->Draw(arcball, projection); }
};
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

    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;
    void textInputEvent(TextInputEvent& event) override;
    void keyPressEvent(KeyEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;
    void updateColor();
    void POM();
    void initGrid();

    template <typename T>
    void dropDownMenu(const char* menuName, T &currentSelectable, std::map<T, std::string> & map);

    void save();
    bool open();
    void saveAs();


    // Tested geometries
    LC::Torus _sheet;
    CrossX _dsheet;

    // Torus with PhongGL shader
    LC::NormalTorus _sheetNormal;

    /*
        GUI Widget
        *** Make sure to keep decoupled from actual simulation core functionality
    */
    Widget _widget;

    LC::Header _header;
    LC::Imaging::UniformGrid::POM _pomImager;
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

    _solver = std::make_unique<LC::FrankOseen::ElasticOnly::FOFDSolver>();

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

    data->k11 = LC::FrankOseen::ElasticConstants::_5CB(LC::FrankOseen::ElasticConstants::Constant::k11);
    data->k22 = LC::FrankOseen::ElasticConstants::_5CB(LC::FrankOseen::ElasticConstants::Constant::k22);
    data->k33 = LC::FrankOseen::ElasticConstants::_5CB(LC::FrankOseen::ElasticConstants::Constant::k33);

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


    initGrid();
    // Colors
    updateColor();


    LC_INFO("Created Sandbox!");
}

Sandbox::~Sandbox() {
	LC_INFO("Destroying Sandbox.");
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

        bool updateImageFromLoad = false;
        bool updateImageFromRadio = false;

        if (_widget.commands.isPressed(KeyEvent::Key::S)) {
            save();
        }

        if (_widget.commands.isPressed(KeyEvent::Key::O)) {
            updateImageFromLoad = open();
        }

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("File")) {

                //if (ImGui::MenuItem("New", "Ctrl+N")) {}
                if (ImGui::MenuItem("Open", "Ctrl+O")) {
                    updateImageFromLoad = open();
                }
                if (ImGui::MenuItem("Save", "Ctrl+S")) {
                    save();
                }
                if (ImGui::MenuItem("Save As..")) {
                    saveAs();
                }

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
            std::map<LC::FrankOseen::LC_TYPE, std::string> lcMap = LC::FrankOseen::LiquidCrystal::Map();
            dropDownMenu<LC::FrankOseen::LC_TYPE>("LC Type", _widget.lcType, lcMap);
        }

        // Cross section radio buttons
        {
            int axesInt = static_cast<int>(_widget.axis);
            
            ImGui::RadioButton("xy", &axesInt, static_cast<int>(Axis::z));
            ImGui::SameLine();
            ImGui::RadioButton("yz", &axesInt, static_cast<int>(Axis::x));
            ImGui::SameLine();
            ImGui::RadioButton("xz", &axesInt, static_cast<int>(Axis::y));

            Axis newAxis = static_cast<Axis>(axesInt);

            if (axesInt != _widget.axis) {
                updateImageFromRadio = true;
                _widget.axis = axesInt;
                initGrid();
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
                    using Waveplate = LC::Imaging::UniformGrid::POM::Waveplate;
                    std::map<Waveplate, std::string> waveplateMap{ { Waveplate::None, "None" }, { Waveplate::Full530nm, "530 nm Full" } };
                    dropDownMenu<Waveplate>("Waveplate", _pomImager.waveplate, waveplateMap);
                }


                // Show additional options in a dropdown (TODO)
                ImGui::InputFloat3("Lamp intensity", &_pomImager.intensity[0]);
                ImGui::InputFloat3("RGB", &_pomImager.lightRGB[0]);
                ImGui::InputFloat("Gamma", &_pomImager.gamma);

                ImGui::Checkbox("Crossed polarizer", &_pomImager.polarizers);
                
                ImGui::InputFloat("Polarizer angle", &_pomImager.polarizerAngle);
            }
        }

        _widget.updateImage = ImGui::Button("Update Image") || updateImageFromLoad
                                                            || updateImageFromRadio;


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


    //_sheet.Draw(_arcballCamera, _projectionMatrix);
    _dsheet.Draw(_arcballCamera, _projectionMatrix);

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


    if (moving || _ioUpdate) redraw();
}

// Todo: inherit all these events through Application
void Sandbox::mousePressEvent(MouseEvent& event) {

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

    // Check if Ctrl + S or Ctrl + O is pressed

    if ((event.key() == KeyEvent::Key::S)) _widget.commands.press(KeyEvent::Key::S);
    else if ((event.key() == KeyEvent::Key::O)) _widget.commands.press(KeyEvent::Key::O);
    else if ((event.key() == KeyEvent::Key::LeftCtrl)) { _widget.commands.ctrl = true; }

    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::keyReleaseEvent(KeyEvent& event) {

    if ((event.key() == KeyEvent::Key::S)) _widget.commands.release(KeyEvent::Key::S);
    else if ((event.key() == KeyEvent::Key::O)) _widget.commands.release(KeyEvent::Key::O);
    else if ((event.key() == KeyEvent::Key::LeftCtrl)) { _widget.commands.ctrl = false; }


    if (_imgui.handleKeyReleaseEvent(event)) _ioUpdate = true;
}

void Sandbox::updateColor() {

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;
    using T4 = LC::FrankOseen::ElasticOnly::FOFDSolver::Tensor4;
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    // Colors
    std::size_t slice = data->voxels[0];
    std::size_t cross_slice = data->voxels[1] * slice;
    std::size_t volslice = data->voxels[2] * cross_slice;

    T4 nn(data->directors, data->voxels[0], data->voxels[1], data->voxels[2], 3);

    Axis ax = static_cast<Axis>(_widget.axis);


    // x -> (1, 2)
    // y -> (0, 2)
    // z -> (0, 1)
    int xx = (ax == Axis::x) ? 1 : 0;
    int yy = (ax == Axis::y) ? 2 : xx + 1;

    std::size_t permutedSlice = data->voxels[xx];

    auto cross_idx = [permutedSlice](int i, int j) {
        return permutedSlice * j + i;
    };

    int hvox = data->voxels[_widget.axis] / 2;

    LC::scalar theta, phi, nx, ny, nz;
    for (int i = 0; i < data->voxels[xx]; i++) {
        for (int j = 0; j < data->voxels[yy]; j++) {

            if (ax == Axis::z) {
                nx = nn(i, j, hvox, 0);
                ny = nn(i, j, hvox, 1);
                nz = nn(i, j, hvox, 2);
            }
            else if (ax == Axis::y) {
                nx = nn(i, hvox, j, 0);
                ny = nn(i, hvox, j, 1);
                nz = nn(i, hvox, j, 2);
            }
            else if (ax == Axis::x) {
                nx = nn(hvox, i, j, 0);
                ny = nn(hvox, i, j, 1);
                nz = nn(hvox, i, j, 2);
            }
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

            _dsheet.section->data[cross_idx(i, j)].color = Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] });
        }
    }

    // Update sheet
    _dsheet.section->sheetBuffer.setData(_dsheet.section->data, GL::BufferUsage::DynamicDraw);
}

// Only for 5CB.
void Sandbox::POM() {

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;
    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    LC::SIscalar pitch = _widget.pitch;
    LC::scalar dop = data->cell_dims[2];

    _pomImager.n0 = LC::FrankOseen::OpticalConstants::LC(_widget.lcType, LC::FrankOseen::OpticalConstants::Constant::n_o).first;
    _pomImager.ne = LC::FrankOseen::OpticalConstants::LC(_widget.lcType, LC::FrankOseen::OpticalConstants::Constant::n_e).first;
 
    _pomImager.thickness = pitch.first * dop;
    _pomImager.dz = _pomImager.thickness * 1e-6 / data->voxels[2]; // meters

    _pomImager.Compute(data->directors, data->voxels, (void*)(&_dsheet.section->data), [](void *data, const std::array<float, 4>&color, std::size_t idx) {
        Magnum::Containers::Array<LC::DynamicColorSheet::Vertex>* dataPtr = (Magnum::Containers::Array<LC::DynamicColorSheet::Vertex>*)data;
        for (int i = 0; i < 4; i++)
            (*dataPtr)[idx].color[i] = color[i];
    });

    // Update sheet
    _dsheet.section->sheetBuffer.setData(_dsheet.section->data, GL::BufferUsage::DynamicDraw);
}

void Sandbox::save() {
    // Check if file exists (was loaded)

    std::string res;

    if (!_widget.loadedFromFile) {
        // File save
        auto sf = pfd::save_file("Select save file name", LCLAB2_ROOT_PATH,
            { "LMT Files (.lmt .lmat)", "*.lmt *.lmat",
              "All Files", "*" },
            pfd::opt::none);


        // Pass file to header to be written
        res = sf.result();
    }
    else {
        res = _header.readFile;
    }

    if (!res.empty()) {

        using FOFDSolver = LC::FrankOseen::ElasticOnly::FOFDSolver;
        FOFDSolver* solver = (FOFDSolver*)_solver.get();

        solver->ConfigureHeader(_header);

        _header.write(res);
        _header.writeBody();
        LC_INFO("Saved to file {0}", res.c_str());
    }
}

bool Sandbox::open() {
    using FOFDSolver = LC::FrankOseen::ElasticOnly::FOFDSolver;
    bool updateImageFromLoad = false;
    
    
    auto of = pfd::open_file("Select save file name", LCLAB2_ROOT_PATH,
                { "LMT Files (.lmt .lmat)", "*.lmt *.lmat",
                        "All Files", "*" },
                pfd::opt::none);

    auto res = of.result();

    if (!res.empty()) {
        FOFDSolver* solver = (FOFDSolver*)_solver.get();

        _header.readFile = res[0];

        solver->ReadDataFromHeader(_header);

        // Reinit grid
        initGrid();

        updateImageFromLoad = true;

        // Tells program to save to readFile if "Save" is pressed
        _widget.loadedFromFile = true;
        LC_INFO("Loaded file {0}", res[0].c_str());

    }

    return updateImageFromLoad;
}

void Sandbox::saveAs() {

    // File save
    auto sf = pfd::save_file("Select save file name", LCLAB2_ROOT_PATH,
        { "LMT Files (.lmt .lmat)", "*.lmt *.lmat",
              "All Files", "*" },
        pfd::opt::none);


    // Pass file to header to be written
    auto res = sf.result();

    if (!res.empty()) {

        using FOFDSolver = LC::FrankOseen::ElasticOnly::FOFDSolver;
        FOFDSolver* solver = (FOFDSolver*)_solver.get();

        solver->ConfigureHeader(_header);

        _header.write(res);
        _header.writeBody();
        LC_INFO("Saved to file {0}", res.c_str());
    }
}

void Sandbox::initGrid() {
    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::Dataset;

    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    int d = _widget.axis;
    Axis ax = static_cast<Axis>(d);
    
    _dsheet.section.reset(new LC::DynamicColorSheet);

    int i = (ax == Axis::x) ? 1 : 0;
    int j = (ax == Axis::y) ? 2 : i + 1;

    _dsheet.section->NX = data->voxels[i];
    _dsheet.section->NY = data->voxels[j];
    _dsheet.section->CX = data->cell_dims[i];
    _dsheet.section->CY = data->cell_dims[j];

    _dsheet.section->Init();
    
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

LC::Application* LC::createApplication(int argc, char **argv) {

	return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}