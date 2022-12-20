#include "DropDownMenu.h"
#include "Widget.h"
#include "Boxes.h"
#include "ZProfile.h"
#include "Quadrilateral.h"
#include "Shader.h"

using namespace Magnum;
using namespace Math::Literals;
using DirectorSolver = LC::FrankOseen::Electric::FOFDSolver;
using AppSolver = LC::Film::Simulation;
using Dataset = DirectorSolver::Dataset;
using Geometry = LC::EllipsoidArray;

struct Surface {

    void Draw(const Matrix4& viewMatrix, const Matrix4& projectionMatrix) {
        Matrix4 mvp = projectionMatrix * viewMatrix;
        for (auto& q : quads)
            q.Draw(mvp);
    }

    void Append(const std::array<Vector3, 4>& verts, const std::array<int, 4>& ids) {
        quads.emplace_back(verts, ids, &shader.program);
    }

    void Append(const std::array<Vector3, 4>& verts) {
        quads.emplace_back(verts, &shader.program);
    }

    void Reserve(unsigned int n) { quads.reserve(n); }

    void Clear() {
        quads.clear();
    }

    Shader shader;
    std::vector<Quadrilateral> quads;
};


#define FULL_BUILD 1

class Sandbox : public LC::Application
{
public:
    explicit Sandbox(const Arguments& arguments);
    ~Sandbox();

    void keyPressEvent(KeyEvent& event) override;

    void updateBoxes(bool first = false);
    void updateFilm();
    void updateFilmGraphics();

    LC::scalar computeStrain();

    void drawBoxes();
    void drawFilm();
    void drawProfile();
    void mainWindow();
    void updateCamera();
    void plotFilmData();
    void initGraphics();

private:
    virtual void drawEvent() override;

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader = Shaders::VertexColorGL3D{};
    SceneGraph::DrawableGroup3D _transparentDrawables;

    Boxes _boxes;
    Surface _filmSurface;
    Widget _widget;

    AppSolver* _sim;

    Containers::Optional<Geometry> _film;
    Containers::Optional<Zprofile> _zprofile;
    std::vector<unsigned int> _surfaceNodes;
    unsigned int _numSurfaceNodes;

    LC::Imaging::UniformGrid::POM _pomImager;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("Film Response Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable) } {
    Utility::Arguments args;
    args.addOption('f', "file", "D:\\dev\\lclab2\\data\\simulations\\homeotropic\\toron_npp30.lmt")
        .setHelp("file", "film file", "F")
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

    // Actual film solver
    _solver = std::make_unique<AppSolver>();
    // Point sim to solver
    _sim = (AppSolver*)(_solver.get());

    // Will be useful when initializing data from a file
    std::unique_ptr<DirectorSolver> directorSolver = std::make_unique<DirectorSolver>();


    // Load a director file
    _header.readFile = "D:\\dev\\lclab2\\data\\simulations\\homeotropic\\smallest_toron.lmt";
    directorSolver->Import(_header);

    LC_INFO("Imported file {0}", _header.readFile.c_str());

    // Make empty again
    {
        //_header.clean();
        //_header.readFile = "";
        //_header.writeFile = "";
        _header = LC::Header{};
    }

    // Populate sim with data
    Dataset* data = (Dataset*)(directorSolver->GetDataPtr());
    // Cell dims stored as d/p
    auto celldims = data->cell_dims;
    LC::scalar cvol = 1.;
    LC::scalar pitch = _widget.pitch * 1e-6; // m

    // cvol = dx * dy * dz
    for (int d = 0; d < 3; d++) {
        celldims[d] *= pitch;
        cvol *= celldims[d] / (data->voxels[d] - 1);
    }


    // Update element masses
    // Note a 10 um x 10 um x 10 um element (vol = 10^-15)
    // has a mass of 10^-12 kg (ref: Joana Aizenberg's paper)
    constexpr LC::scalar base_vol = 1e-15;

    // Initialize the film graphics
    _film = Geometry{};

#if FULL_BUILD
    _sim->Parameters().bulkMass *= cvol / base_vol;
    _sim->Parameters().surfaceMass *= cvol / base_vol;
    _sim->Parameters().edgeMass *= cvol / base_vol;
    _sim->Parameters().cornerMass *= cvol / base_vol;
    _sim->Parameters().PBCs = true;
    _sim->InitVerticesAndElements(data->directors.get(), data->voxels, celldims);

#else // testing
    {
        _sim->Parameters().PBCs = false;
        _sim->Parameters().dT = 10e-9; // 10 ns
        LC::scalar UdS = -570; // kPa
        _sim->Parameters().UdS = UdS * 1e3; // Pa
        _sim->Parameters().directorAxis = { 0., 0., 1. };
        std::array<int, 3> vox = { 2,2,2 };
        _sim->InitVerticesAndElements(0, vox, { 1e-7, 1e-7, 1e-7 });
        _film->scale = 0.01;
        _widget.drawZProfile = false;
        _widget.drawFilm = true;
        _widget.globalStrain = true;
        _widget.strain_recorder = StrainRecorder(UdS, UdS, 50);
        _widget.strain_recorder.tolerance = 1e-8;
        _widget.strain_recorder.film_thickness = _sim->Parameters().cell_dims[2] * 1e6;
    }
#endif

    initGraphics();

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


void Sandbox::updateCamera() {

    /* Update camera */
    if (_cameraType == CameraType::ArcBall) {
        _arcballCamera->updateTransformation();

        polyRenderer();
        updateBoxes();
        drawBoxes();

        if (_widget.updateFilm)
            updateFilm();

        if (_widget.drawFilm)
            drawFilm();

        if (_widget.drawZProfile)
            drawProfile();

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

            initGraphics();
        }
    }
    if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
}

void Sandbox::initGraphics() {
    // Function to access data
    std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* data, std::size_t i) {
        Magnum::Vector3* vdata = (Magnum::Vector3*)data;
        return vdata[i];
    };

    LC::Film::Vertex* nodes = _sim->getNodes();
    _numSurfaceNodes = 0;

    LC::Film::Vertex::Role role = LC::Film::Vertex::Role::Surface;

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
            // Convert to um
            tempNodes[count] = 1e6f * Vector3(nodes[i].position[0], nodes[i].position[1], nodes[i].position[2]);

            ++count;
        }
    }

    _film->Init(&tempNodes[0], access, _numSurfaceNodes);

    // Populate the z profile
    _zprofile = Zprofile{};

    updateFilmGraphics();

    updateBoxes();
}

void Sandbox::updateBoxes(bool first) {

    _boxes.Clear();

    auto dims = _sim->Parameters().cell_dims;
    Vector3 s;
    
    s[0] = dims[0] * 0.5e6f;
    s[1] = dims[1] * 0.5e6f;
    s[2] = dims[2] * 0.5e6f;


    // Draw a 1x1x1 cube to the screen for reference
    _boxes.Append(_arcballCamera->viewMatrix() *
         Matrix4::scaling(s), 0xff0000_rgbf);
}

LC::scalar Sandbox::computeStrain() {
    
    int vox_x = _sim->Parameters().voxels[0];
    int vox_y = _sim->Parameters().voxels[1];
    int vox_z = _sim->Parameters().voxels[2];
    int slice = vox_x * vox_y;
    auto index = [&](int x, int y, int z) {
        return x + vox_x * y + slice * z;
    };

    LC::scalar zmax = _sim->getNodes()[index(0, 0, vox_z - 1)].position.z();
    LC::scalar zmin = zmax;
    for (int i = 0; i < vox_x; i++) {
        for (int j = 0; j < vox_y; j++) {
            
            LC::scalar znew = _sim->getNodes()[index(i, j, vox_z - 1)].position.z();
            if (zmax < znew) zmax = znew;
            if (zmin > znew) zmin = znew;
        }
    }
    
    LC::scalar strain;
    if (_widget.globalStrain)
        strain = 100. * (zmax - 0.5 * _sim->Parameters().cell_dims[2]) / _sim->Parameters().cell_dims[2];
    else
        strain = 100. * (zmax - zmin) / _sim->Parameters().cell_dims[2];
    return strain;
}

void Sandbox::updateFilm() {

#if FULL_BUILD
    for (int i = 0; i < _widget.updateCycle; i++)
        _sim->updateNodes();
#else
    // Update until strain stops changing
    // Compute strain
    LC::scalar strain_current = computeStrain();
    LC::scalar strain_last = strain_current - 1.;
    while (abs(strain_last - strain_current) > _widget.strain_recorder.tolerance) {
        _sim->updateNodes();
        strain_last = strain_current;
        strain_current = computeStrain();
    }

#endif

    if (_widget.updateGraphics) {
        updateFilmGraphics();
    }
    // Film has been updated
#if FULL_BUILD
    if (!_widget.continuousUpdate)
        _widget.updateFilm = false;
#else

    if (_widget.strain_recorder.Finished()) {
        _widget.updateFilm = false;
        // Write data
        _widget.strain_recorder.Write();
    }
    else {
        _widget.strain_recorder.Add(strain_current);
        // Update UdS if not done
        if (!_widget.strain_recorder.Finished())
            _sim->Parameters().UdS += _widget.strain_recorder.change_in_UdS * 1e3;
    }

#endif
}


void Sandbox::updateFilmGraphics() {

    Float polyRadius = _film->scale / (Float)pow(_film->numObjects, 1.0f / 3.0f);

    Matrix4 rotation;
    Float hPi = M_PI / 2.0;

    auto nodes = _sim->getNodes();
    auto voxels = _sim->Parameters().voxels;

    for (int i = 0; i < _numSurfaceNodes; i++) {

        float theta = acos(nodes[_surfaceNodes[i]].director[2]);
        float phi = atan2(nodes[_surfaceNodes[i]].director[1], nodes[_surfaceNodes[i]].director[0]);

        _film->polyPositions[i] = 1e6 * Vector3(nodes[_surfaceNodes[i]].position[0],
            nodes[_surfaceNodes[i]].position[1],
            nodes[_surfaceNodes[i]].position[2]);

        rotation = Matrix4::rotationZ(Rad{ -hPi + phi }) * Matrix4::rotationX(Rad{ hPi - theta });
        _film->polyInstanceData[i].transformationMatrix =
            Matrix4::translation(_film->polyPositions[i]) * Matrix4::scaling(Vector3{ polyRadius }) * rotation;
        _film->polyInstanceData[i].normalMatrix =
            _film->polyInstanceData[i].transformationMatrix.normalMatrix();
    }

    // Z profile
    std::unique_ptr<LC::scalar[]> rz(new LC::scalar[voxels[0] * voxels[1]]);
    for (int i = 0; i < voxels[0]; i++)
        for (int j = 0; j < voxels[1]; j++) {
            unsigned int idx = i + voxels[0] * j;
            rz[idx] = (_sim->getNodes() + (voxels[2] - 1) * voxels[0] * voxels[1] + idx)->position(2) * 1e6;
        }

    // Convert to um
    auto cdims = _sim->Parameters().cell_dims;
    for (auto &c : cdims)
        c *= 1e6;
    
    // Position rz has units um like cell_dims
    _zprofile->GenerateProfile(rz.get(), voxels, cdims);

    // Create the film surface

    _filmSurface.Clear();


    auto voxElems = voxels;
    for (auto& v : voxElems)
        v -= 1;

    unsigned int nQuads = 2 * (voxElems[0] * voxElems[1] + voxElems[0] * voxElems[2] + voxElems[2] * voxElems[1]);

    _filmSurface.Reserve(nQuads);

    std::array<Vector3, 4> verts;


    // +-x faces
    for (int i : { 0, voxElems[0] }) {

        if (i == 0) {
            if (!_widget.surface_selector.axes[0].draw[0])
                continue;
        }
        else if (!_widget.surface_selector.axes[0].draw[1])
            continue;

        for (int j = 0; j < voxElems[1]; j++)
            for (int k = 0; k < voxElems[2]; k++) {


                auto v1 = 1e6 * nodes[i + j * voxels[0] + k * voxels[0] * voxels[1]].position;
                auto v2 = 1e6 * nodes[i + (j + 1) * voxels[0] + k * voxels[0] * voxels[1]].position;
                auto v3 = 1e6 * nodes[i + (j + 1) * voxels[0] + (k + 1) * voxels[0] * voxels[1]].position;
                auto v4 = 1e6 * nodes[i + j * voxels[0] + (k + 1) * voxels[0] * voxels[1]].position;

                verts[0] = Vector3(v1(0), v1(1), v1(2));
                verts[1] = Vector3(v2(0), v2(1), v2(2));
                verts[2] = Vector3(v3(0), v3(1), v3(2));
                verts[3] = Vector3(v4(0), v4(1), v4(2));

                _filmSurface.Append(verts);
            }
    }

    // +-y faces
    for (int j : { 0, voxElems[0] }) {

        if (j == 0) {
            if (!_widget.surface_selector.axes[1].draw[0])
                continue;
        }
        else if (!_widget.surface_selector.axes[1].draw[1])
            continue;

        for (int i = 0; i < voxElems[0]; i++)
            for (int k = 0; k < voxElems[2]; k++) {




                auto v1 = 1e6 * nodes[i + j * voxels[0] + k * voxels[0] * voxels[1]].position;
                auto v2 = 1e6 * nodes[(i + 1) + j * voxels[0] + k * voxels[0] * voxels[1]].position;
                auto v3 = 1e6 * nodes[(i + 1) + j * voxels[0] + (k + 1) * voxels[0] * voxels[1]].position;
                auto v4 = 1e6 * nodes[i + j * voxels[0] + (k + 1) * voxels[0] * voxels[1]].position;

                verts[0] = Vector3(v1(0), v1(1), v1(2));
                verts[1] = Vector3(v2(0), v2(1), v2(2));
                verts[2] = Vector3(v3(0), v3(1), v3(2));
                verts[3] = Vector3(v4(0), v4(1), v4(2));

                _filmSurface.Append(verts);
            }
    }

    // +-z faces
    for (int k : { 0, voxElems[2] }) {

        if (k == 0) {
            if (!_widget.surface_selector.axes[2].draw[0])
                continue;
        }
        else if (!_widget.surface_selector.axes[2].draw[1])
            continue;

        for (int i = 0; i < voxElems[0]; i++)
            for (int j = 0; j < voxElems[1]; j++) {

                auto v1 = 1e6 * nodes[i + j * voxels[0] + k * voxels[0] * voxels[1]].position;
                auto v2 = 1e6 * nodes[i + 1 + j * voxels[0] + k * voxels[0] * voxels[1]].position;
                auto v3 = 1e6 * nodes[i + 1 + (j + 1) * voxels[0] + k * voxels[0] * voxels[1]].position;
                auto v4 = 1e6 * nodes[i + (j + 1) * voxels[0] + k * voxels[0] * voxels[1]].position;

                verts[0] = Vector3(v1(0), v1(1), v1(2));
                verts[1] = Vector3(v2(0), v2(1), v2(2));
                verts[2] = Vector3(v3(0), v3(1), v3(2));
                verts[3] = Vector3(v4(0), v4(1), v4(2));

                _filmSurface.Append(verts);
            }
    }

}

void Sandbox::drawFilm() {
    _film->Draw(_arcballCamera, _projectionMatrix);
    _filmSurface.Draw(_arcballCamera->viewMatrix(), _projectionMatrix);
}

void Sandbox::drawProfile() {
    if (_zprofile)
        _zprofile->Draw(_arcballCamera, _projectionMatrix);
}

void Sandbox::drawBoxes() {
    _boxes.DrawFrame(_projectionMatrix);
}

void Sandbox::mainWindow() {

    if (!_widget.showSettings) return;

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar;
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::Begin("film activation", 0, window_flags);

    // Simulation details
    {
        std::string pitch = "pitch = " + std::to_string(_widget.pitch) + " um";
        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, pitch.c_str());
        std::string cdims = "Cell dimensions = ("
            + std::to_string(_sim->Parameters().cell_dims[0] * 1e6)
            + " um X "
            + std::to_string(_sim->Parameters().cell_dims[1] * 1e6)
            + " um X "
            + std::to_string(_sim->Parameters().cell_dims[2] * 1e6)
            + " um)";
        std::string vox = "Voxels = ("
            + std::to_string(_sim->Parameters().voxels[0])
            + " X "
            + std::to_string(_sim->Parameters().voxels[1])
            + " X "
            + std::to_string(_sim->Parameters().voxels[2])
            + ")";
        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, cdims.c_str());
        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, vox.c_str());


        std::string PBC = "PBCs: ";
        if (_sim->Parameters().PBCs) PBC += "on";
        else PBC += "off";

        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, PBC.c_str());
        std::string iterations = "Iterations relaxed = " + std::to_string(_sim->Parameters().iterations);
        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, iterations.c_str());

        std::string strain = "Strain = " + std::to_string(computeStrain()) + "%%";
        ImGui::TextColored({ 0.f, 1.f, 0.f, 1.f }, strain.c_str());
    }

    ImGui::SliderFloat("Node scale", &_film->scale, 0.125f, 10.0f);
    ImGui::SliderInt("Chain extensibility", &_sim->Parameters().I_m, 4, 30);

    float kappa = _sim->Parameters().kappa * 1e-3;
    ImGui::SliderFloat("Bulk modulus (kPa)", &kappa, 0.f, 1e5f);
    _sim->Parameters().kappa = kappa * 1e3;

    float mu = _sim->Parameters().mu * 1e-3;
    ImGui::SliderFloat("Shear modulus (kPa)", &mu, 0.f, 800.f);
    _sim->Parameters().mu = mu * 1e3;

    float UdS = _sim->Parameters().UdS *1e-3;
    ImGui::SliderFloat("UdS (kPa)", &UdS, -800.f, 800.f);
    _sim->Parameters().UdS = UdS * 1e3;

    float dT = _sim->Parameters().dT * 1e9;
    ImGui::InputFloat("Time step (ns)", &dT);
    _sim->Parameters().dT = dT * 1e-9;

    ImGui::InputInt("Update iterations", &_widget.updateCycle);
    

    if (ImGui::Button("Update film")) {
        _widget.updateFilm = true;

    }

    ImGui::SameLine();

    ImGui::Checkbox("Continuous update", &_widget.continuousUpdate);

    ImGui::SameLine();

    ImGui::Checkbox("Global strain", &_widget.globalStrain);

    ImGui::TextColored({ 1.f, 1.f, 0.f, 1.f }, "Draw");
    ImGui::Checkbox("Draw film", &_widget.drawFilm);
    ImGui::SameLine();
    ImGui::Checkbox("Draw z-profile", &_widget.drawZProfile);
    ImGui::Checkbox("Update graphics", &_widget.updateGraphics);

    _widget.surface_selector.DrawGUI();

    // Show the film
    plotFilmData();


    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
        1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
    ImGui::End();
}

void Sandbox::plotFilmData() {
    auto voxels = _sim->Parameters().voxels;

    if (ImGui::RadioButton("x-line", &_widget.radioZProfileAxis, 0))
        _widget.radioZProfileIndex = -1;
    ImGui::SameLine();
    if (ImGui::RadioButton("y-line", &_widget.radioZProfileAxis, 1))
        _widget.radioZProfileIndex = -1;
    ImGui::SameLine();
    ImGui::PushItemWidth(150.0f);
    if (_widget.radioZProfileIndex == -1)
        _widget.radioZProfileIndex = (voxels[_widget.radioZProfileAxis] - 1) / 2;
    ImGui::SliderInt("Index", &_widget.radioZProfileIndex, 0, voxels[_widget.radioZProfileAxis] - 1);
    ImGui::PopItemWidth();

    if (_zprofile) {
        // Plot of the z-profile along x slice

        int vox[] = { _zprofile->graph->vox_x, _zprofile->graph->vox_x };
        int complementAxis = (_widget.radioZProfileAxis + 1) % 2;

        int vox_x = vox[_widget.radioZProfileAxis];
        int vox_y = vox[complementAxis];

        std::unique_ptr<float[]> ax(new float[vox_x]);
        std::unique_ptr<float[]> line_scan(new float[vox_y]);

        // Convert to um
        auto cdims = _sim->Parameters().cell_dims;
        for (auto& c : cdims)
            c *= 1e6;

        auto voxels = _sim->Parameters().voxels;

        float dx = cdims[_widget.radioZProfileAxis] / (voxels[_widget.radioZProfileAxis] - 1);
        float dy = cdims[complementAxis] / (voxels[complementAxis] - 1);

        for (int i = 0; i < vox_x; i++) {
            if (_widget.radioZProfileAxis == 0) {
                ax[i] = _zprofile->graph->data[i + _widget.radioZProfileIndex * vox[0]].x(); // um
                line_scan[i] = _zprofile->graph->data[i + _widget.radioZProfileIndex * vox[0]].z() * 1e3; // nm
            }
            else if (_widget.radioZProfileAxis == 1) {
                ax[i] = _zprofile->graph->data[_widget.radioZProfileIndex + i * vox[0]].y(); // um
                line_scan[i] = _zprofile->graph->data[_widget.radioZProfileIndex + i * vox[0]].z() * 1e3; // nm
            }
        }

        std::string letter;

        if (_widget.radioZProfileAxis == 0) // x-axis
            letter = 'x';
        else
            letter = 'y';

        std::string plt_title = "Z-profile plot (" + letter + "-axis)";

        auto fautofit = ImPlotAxisFlags_::ImPlotAxisFlags_AutoFit;

        if (ImPlot::BeginPlot(plt_title.c_str(), "L (um)", "h (nm)", ImVec2(-1, 0), 0, fautofit, fautofit)) {
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
            ImPlot::PlotLine(("Z(" + letter + ")").c_str(), ax.get(), line_scan.get(), vox_x);
            ImPlot::EndPlot();
        }

        // Draw the plane corresponding to the cross section

        _zprofile->box.dims[complementAxis] = 0;
        _zprofile->box.dims[_widget.radioZProfileAxis] = cdims[_widget.radioZProfileAxis];
        _zprofile->box.dims[2] = cdims[2];

        _zprofile->box.translation[complementAxis] = dy * (_widget.radioZProfileIndex - voxels[complementAxis] + 1) + 0.5f * cdims[complementAxis];
        _zprofile->box.translation[_widget.radioZProfileAxis] = 0.0f;
        _zprofile->box.translation[2] = 0.0f;

    }
}

LC::Application* LC::createApplication(int argc, char** argv) {

    return new Sandbox{ Platform::Application::Arguments{argc, argv} };
}