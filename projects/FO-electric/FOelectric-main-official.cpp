#include "Widget.h"
#include "find_components.h"
#include "VKHelpers.h"
#include "ZProfile.h"
#include "LehmanCluster.h"
#include "CrossSection.h"
#include "SmoothIsosurface.h"
#include "QuaternionFromBasis.h"
#include "NematicMultiplaneManager.h"
#include <tuple>
#include <filesystem>
#include <algorithm>
#include <set>

using FOFDSolver = LC::FrankOseen::Electric::FOFDSolver;
using Dataset = FOFDSolver::Dataset;

// Export vortex knots for KnotPlot
std::tuple<std::vector<Eigen::Vector3d>, std::vector<int>> exportVortexKnot(const std::vector<MeshLib::PNCVertex<float>>& vertices,
    const std::vector<MeshLib::Triangle>& triangles,
    LC::scalar min_point_dist, const std::string& fname, bool append = false);

// Determine heliknoton COM by preimage
std::vector<Eigen::Vector3d> detect_heliknoton_COM(std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>>&,const std::vector<Eigen::Vector3d>&);

struct BoxInstanceData {
    Matrix4 transformationMatrix;
    Color3 color;
};

struct SphereInstanceData {
    Matrix4 transformationMatrix;
    Color3 color;
};

// Base properties to inherit for a preimage
struct PreimageBase {
    bool draw = true;
    float isoLevel = 0.07f;
    float alpha = 1.0f;
    bool opentab = true;
    LC::Surface surface;

    Containers::Optional<GL::Mesh> mesh;
};

struct S2Fiber : public PreimageBase {
    S2Fiber() {}
    S2Fiber(float t, float p, float iso = 0.07f) : theta(t), phi(p) { isoLevel = iso; }
    
    friend bool operator == (const S2Fiber& p1, const S2Fiber& p2) {
        return (p1.theta == p2.theta && p1.phi == p2.phi) ? 1 : 0;
    }

    float theta = 0.0f;
    float phi = 0.0f;
    int subdivisions = 0;
    bool normal_inversion = false;
    bool domain_inversion = false;
    bool cullFaces = true;
    bool pontryagin = false;
    bool draw = true;
};

struct PionPreimage : PreimageBase {
    PionPreimage() = default;
    PionPreimage(float p1, float p2, float p3, float iso = 0.07f) : pi1(p1), pi2(p2), pi3(p3) { isoLevel = iso; }
    PionPreimage(const std::array<float, 3> &pi, float iso = 0.07f) : pi1(pi[0]), pi2(pi[1]), pi3(pi[2]) { isoLevel = iso; }

    friend bool operator == (const PionPreimage& p1, const PionPreimage& p2) {
        return (p1.pi1 == p2.pi1 && p1.pi2 == p2.pi2 && p1.pi3 == p2.pi3) ? 1 : 0;
    }

    float pi1 = 0.0f;
    float pi2 = 0.0f;
    float pi3 = 1.0f;
    bool noCull = false;
};

struct VortexLine {
    LC::Math::Interpolant_d_1D interpolant;
    std::unique_ptr<std::size_t[]> neighbors;
    std::unique_ptr<LC::scalar[]> positions;
    std::unique_ptr<std::size_t[]> queryDomain;
    std::size_t numNodes;
    std::size_t knn = 35;
    int maxPoints = 150;
    int maxConic = 30;
    // 3 times the background grid spacing
    float pointSeparation = 1.5f;
    // Number of times line is post processed with Chaikin's method
    int upSample = 2;
    float tubeRadius = 1.f;

    // Generated data
    std::vector<Eigen::Vector3d> components_COM;
    std::vector<std::vector<uint>> components;
    std::vector<MeshLib::PNCVertex<float>> vertices;
    std::unique_ptr<short[]> valid_field;
};

struct VortexKnot : public PreimageBase {
    VortexKnot(float iso = 0.07f) { isoLevel = iso; }

    friend bool operator == (const VortexKnot& v1, const VortexKnot& v2) {
        return (v1.isoLevel == v2.isoLevel) ? 1 : 0;
    }

    void UpdateColor() {
        for (auto& v : surface.vertices) {
           v.color = knotColor;
        }
        
        mesh = surface.Mesh();
    }

    Color4 knotColor{ 1.f, 0.f, 0.f, 1.f };
    bool cullFaces = false;
    bool normalInversion = true;
    bool ribbon = false;
    LC::SphereArray points;
};

struct VortexShell : public PreimageBase {
    VortexShell(float iso = 0.07f) { isoLevel = iso; alpha = 0.5f; }

    friend bool operator == (const VortexShell& v1, const VortexShell& v2) {
        return (v1.isoLevel == v2.isoLevel) ? 1 : 0;
    }

    bool invertNormals = false;
    bool noCull = false;
    float queryValue = 0.0f;
};

struct BaryonIsosurface : public PreimageBase {
    BaryonIsosurface() {
        alpha = 0.5f;
        isoLevel = 1.999e-2f;
    }
    friend bool operator == (const BaryonIsosurface& v1, const BaryonIsosurface& v2) {
        return (v1.isoLevel == v2.isoLevel) ? 1 : 0;
    }
    float query = 2e-2f;
    bool noCull = false;
};

struct NematicVisual {
    LC::EllipsoidArray lambda, chi, tau;
    void Draw(const Magnum::Matrix4& viewMatrix, const Magnum::Matrix4& projectionMatrix) {
        if (drawlambda)
            lambda.Draw(viewMatrix, projectionMatrix);
        if (drawchi)
            chi.Draw(viewMatrix, projectionMatrix);
        if (drawtau)
            tau.Draw(viewMatrix, projectionMatrix);
    }
    bool drawlambda = false;
    bool drawchi = false;
    bool drawtau = false;
    int dim1, dim2;
};


class Sandbox : public LC::Application {
public:

    explicit Sandbox(const Arguments& arguments);
    ~Sandbox();

private:
    virtual void drawEvent() override;
    void keyPressEvent(KeyEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;
    void updateColor();
    void POM(std::string override_name = "");
    void initVisuals();
    void updateBoxes(bool first = false);
    void updateSpheres(bool first = false);
    void drawBoxes();
    void drawSpheres();
    void imageMenu();
    void handlePOMWindow();
    void handleVisualPropertiesWindow();
    void handlePreimageWindow();
    void handleVortexKnotWindow();
    void handleKnotInteractionWindow();
    void nematicPreimageManager();
    void handleNonlinearImagingWindow();
    void handleLCINFOWindow();
    void handleModificationWindow();
    void handleInteractionEnergyWindow();
    void handleZProfileWindow();
    void handleMultiplaneWindow();
    void handleLehmannWindow();
    std::tuple <std::vector<MeshLib::PNCVertex<float>>, std::vector<MeshLib::Triangle>>
        TaubinSmoothingGeneral(float lambda,
            float mu,
            int smoothingIterations,
            LC::Math::IsoVertex* verts,
            unsigned int nVert,
            unsigned int* indices,
            unsigned int nInd,
            const float* field_nn,
            const std::array<int, 3>& Vox,
            const std::array<float, 4>& color,
            bool invert_normals,
            int color_style,
            unsigned int nSubdivisions,
            bool saveObj = false,
            const std::string& obj_name = "");
    std::tuple <std::vector<MeshLib::PNCVertex<float>>, std::vector<MeshLib::Triangle>>
        TaubinSmoothing(LC::Math::IsoVertex* verts,
        unsigned int nVert,
        unsigned int* indices,
        unsigned int nInd,
        const float* field_nn,
        const std::array<int, 3>& Vox,
        const std::array<float, 4>& color,
        bool invert_normals,
        int color_style,
        unsigned int nSubdivisions,
        bool saveObj = false,
        const std::string& obj_name = "");
    void generateIsosurface();
    void computeEnergy();
    void repeatVolume(int i);

    void findVortexKnotComponents(bool saveObj = false, const std::string& obj_name = "test.obj");

    std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>> findVortexKnot(bool saveObj, const std::string& obj_name,
        std::unique_ptr<float[]>& chi_field, std::unique_ptr<float[]>& field_nn, std::unique_ptr<float[]>& valid_fieldf,
        std::unique_ptr<float[]>& field_S, std::unique_ptr<float[]>& sample_grid,
        std::unique_ptr<float[]>& field_Q, std::unique_ptr<float[]>& field_Ssb,
        LC::ExtendedMC::MarchingCubes & mc, LC::ExtendedMC::MarchingCubes &mc2,
        std::array<int,3> &vox_interp);

    template <typename T>
    void dropDownMenu(const char* menuName, T& currentSelectable, std::map<T, std::string>& map);

    // Used to draw CrossX mesh
    Shaders::VertexColorGL3D _transparentShader;
    Shaders::PhongGL _phongShader;

    SceneGraph::DrawableGroup3D _transparentDrawables;
    SceneGraph::DrawableGroup3D _transparentNormalDrawables;

    std::unique_ptr<LC::Drawable::Object3D> _preimageManipulator;
    std::unique_ptr<LC::Drawable::Object3D> _vortexManipulator, _processedVortexManipulator;

    Containers::Array<CrossX> _crossSections;

    GL::Mesh _sphereMesh{ NoCreate };
    GL::Buffer _sphereInstanceBuffer{ NoCreate };
    Containers::Array<SphereInstanceData> _sphereInstanceData;
    GL::Mesh _boxMesh{ NoCreate };
    GL::Buffer _boxInstanceBuffer{ NoCreate };
    Shaders::FlatGL3D _boxShader{ NoCreate };
    Containers::Array<BoxInstanceData> _boxInstanceData;
    Containers::Optional<Zprofile> _zprofile;
    std::unique_ptr<LC::scalar[]> _baryon_density;

    /*
        GUI Widget
        *** Make sure to keep decoupled from actual simulation core functionality
    */
    Widget _widget;
    LC::Header _header;
    LC::Imaging::UniformGrid::POM _pomImager;
    LC::Math::Isosurface<float*, float> _isoGenerator;
    std::list<S2Fiber> _nematicPreimages;
    std::list<VortexKnot> _vortexKnot, _processedVortexKnot;
    Containers::Optional<VortexShell> _vortexShell;
    Containers::Optional<PionPreimage> _pionPreimage;
    Containers::Optional<BaryonIsosurface> _baryonIsosurface;
    VortexLine _vortex_line;

    lehman_cross_section _lehman_xsection;

    std::unique_ptr<Eigen::Quaternion<LC::scalar>[]> _quaternion_field;
    Containers::Optional<NematicVisual> _quaternion_plane;
    NematicMultiplaneManager _multiplane_manager;

    LC::Imaging::ImageSeries _image_series, _nonlin_image_series;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("FO Electric Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable),
                                                                GLConfiguration{}.setSampleCount(16)

} {

    Utility::Arguments args;
    args.addOption('q', "topological-charge", "1")
        .setHelp("topological-charge", "topological charge", "Q")
        .addOption('n', "nodes-per-pitch", "15")
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
        .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB());

    // Configuration needs voxels and cell dims
    (*data).Configuration(Dataset::Heliknoton(Q, data->voxels, data->cell_dims, 1., 1.));

    /*
        Set voltage
    */

    (*data).vconfig = [v0](FOFDSolver::Tensor3& vv, int i, int j, int k, int* voxels) {
        vv(i, j, k) = v0 * (k / (voxels[2] - 1));
    };

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
    
    _sphereInstanceBuffer = GL::Buffer{};
    _sphereMesh = MeshTools::compile(Primitives::icosphereWireframe());
    _sphereMesh.addVertexBufferInstanced(_sphereInstanceBuffer, 1, 0,
        Shaders::FlatGL3D::TransformationMatrix{},
        Shaders::FlatGL3D::Color3{});
    
    _widget.heliknotonInsertLocation = { data->voxels[0] / 2, data->voxels[1] / 2, data->voxels[2] / 2 };
    
    updateSpheres(true);

    _widget.series_x_axis_vec.resize(_widget.max_graph_points);
    _widget.energy_series_vec.resize(_widget.max_graph_points);

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
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    _imgui.newFrame();
    Dataset* data = static_cast<Dataset*>(_solver->GetDataPtr());
    FOFDSolver* solver = static_cast<FOFDSolver*>(_solver.get());

    auto normalizeLoadedDirectors = [this]() {
        initVisuals();
        Dataset* dataloc = static_cast<Dataset*>(_solver->GetDataPtr());

        size_t grid_vol = static_cast<size_t>(dataloc->voxels[0]) * dataloc->voxels[1] * dataloc->voxels[2];
        size_t invalid_dirs = 0;

        for (size_t i = 0; i < grid_vol; ++i) {
            LC::scalar nx = dataloc->directors[i];
            LC::scalar ny = dataloc->directors[i + grid_vol];
            LC::scalar nz = dataloc->directors[i + 2 * grid_vol];
            LC::scalar mag = std::sqrt(nx * nx + ny * ny + nz * nz);

            if (mag != 0.) {
                dataloc->directors[i] /= mag;
                dataloc->directors[i + grid_vol] /= mag;
                dataloc->directors[i + 2 * grid_vol] /= mag;
            }
            else {
                ++invalid_dirs;
            }
        }

        if (invalid_dirs > 0) {
            LC_WARN("Loaded degenerate director field");
        }

        _widget.energy_series.clear();
        _widget.series_x_axis.clear();
        for (size_t i = 0; i < _widget.energy_series_vec.size(); ++i) {
            _widget.energy_series_vec[i] = 0.0;
            _widget.series_x_axis_vec[i] = dataloc->numIterations + static_cast<int>(i);
        }

        _image_series = LC::Imaging::ImageSeries(
            static_cast<std::int32_t>(dataloc->voxels[0]),
            static_cast<std::int32_t>(dataloc->voxels[1]),
            _header.readFile + "-POM"
        );
        for (int d = 0; d < 3; ++d) {
            _widget.celldims[d] = dataloc->cell_dims[d];
        }
        };

    if (_widget.showSettings) {
        ImGui::SetNextWindowPos(ImVec2(0.f, 0.f), ImGuiCond_Always);
        ImGui::SetNextWindowSizeConstraints(ImVec2(370.f, 400.f), ImVec2(520.f, ImGui::GetIO().DisplaySize.y));

        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 4.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_GrabRounding, 4.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_PopupRounding, 4.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(6.0f, 4.0f));

        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_AlwaysAutoResize;
        ImGui::Begin("lclab2", nullptr, window_flags);

        saveMenu(_widget.updateImageFromLoad, normalizeLoadedDirectors);
        handleVisualPropertiesWindow();
        imageMenu();

        float avail_w = ImGui::GetContentRegionAvail().x;
        float half_btn_w = (avail_w - ImGui::GetStyle().ItemSpacing.x) * 0.5f;

        if (ImGui::Button("LC Info", ImVec2(half_btn_w, 0.f))) {
            _widget.showLCINFO ^= true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Modify LC", ImVec2(half_btn_w, 0.f))) {
            _widget.showModificationWindow ^= true;
        }

        ImGui::Spacing();

        if (ImGui::CollapsingHeader("Material & Field Parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
            auto lcMap = LC::FrankOseen::LiquidCrystal::Map();
            dropDownMenu<LC::FrankOseen::LC_TYPE>("LC Type", data->lc_type, lcMap);

            if (data->lc_type != LC::FrankOseen::LC_TYPE::CUSTOM) {
                (*data).ElasticConstants(data->lc_type)
                    .ElectricConstants(data->lc_type)
                    .OpticalConstants(data->lc_type);
            }

            bool modified = false;
            ImGui::Spacing();
            ImGui::TextColored(ImVec4(0.4f, 0.85f, 1.0f, 1.0f), "Elastic, Dielectric & Optical Constants");

            if (ImGui::BeginTable("LCParamsTable", 3, ImGuiTableFlags_SizingStretchSame)) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                modified |= ImGui::InputDouble("##k11", &data->k11.first, 0.0, 0.0, "k11: %.3f");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                modified |= ImGui::InputDouble("##k22", &data->k22.first, 0.0, 0.0, "k22: %.3f");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                modified |= ImGui::InputDouble("##k33", &data->k33.first, 0.0, 0.0, "k33: %.3f");

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                modified |= ImGui::InputDouble("##epar", &data->epar, 0.0, 0.0, "epar: %.2f");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                modified |= ImGui::InputDouble("##eper", &data->eper, 0.0, 0.0, "eper: %.2f");
                ImGui::TableNextColumn();
                ImGui::Dummy(ImVec2(0.f, 0.f));

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                modified |= ImGui::InputDouble("##n0", &data->n0, 0.0, 0.0, "n0: %.3f");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                modified |= ImGui::InputDouble("##ne", &data->ne, 0.0, 0.0, "ne: %.3f");
                ImGui::TableNextColumn();
                ImGui::Dummy(ImVec2(0.f, 0.f));

                ImGui::EndTable();
            }

            if (modified) {
                data->lc_type = LC::FrankOseen::LC_TYPE::CUSTOM;
            }

            ImGui::Separator();
            ImGui::PushItemWidth(120.0f);
            ImGui::InputFloat("Voltage (V)", &_widget.voltage, 0.1f, 1.0f, "%.2f");
            ImGui::InputInt("Prerelax Iters", &_widget.voltage_iterations);
            ImGui::PopItemWidth();

            if (ImGui::Button("Apply Voltage", ImVec2(-FLT_MIN, 0.f))) {
                solver->SetVoltage(_widget.voltage, _widget.voltage_iterations);
                LC_INFO("Updated voltage");
            }
            ImGui::Spacing();
        }

        if (ImGui::CollapsingHeader("Cross-Sections & Camera", ImGuiTreeNodeFlags_DefaultOpen)) {
            bool prior_midplane = _widget.midplane;
            ImGui::Checkbox("Lock to Midplane", &_widget.midplane);
            if (!_widget.midplane) {
                if (prior_midplane) {
                    for (int d = 0; d < 3; ++d) {
                        _widget.iPlane[d] = data->voxels[d] / 2;
                    }
                }
                int max_vox = *std::max_element(data->voxels.begin(), data->voxels.end());
                ImGui::SliderInt3("Slice Index", &_widget.iPlane[0], 0, max_vox);
                for (int d = 0; d < 3; ++d) {
                    _widget.iPlane[d] = std::clamp(_widget.iPlane[d], 0, data->voxels[d] - 1);
                }
            }

            if (ImGui::Button("Planes...")) ImGui::OpenPopup("plane_selector");
            ImGui::SameLine();
            if (ImGui::Button("Nematic...")) ImGui::OpenPopup("nematic_plane_selector");
            ImGui::SameLine();
            ImGui::Checkbox("Coupled", &_widget.couplePlaneAndNematic);

            constexpr std::array<const char*, 3> plane_names = { "yz (x-axis)", "xz (y-axis)", "xy (z-axis)" };

            if (ImGui::BeginPopup("plane_selector")) {
                ImGui::TextDisabled("Active Cross-Sections");
                ImGui::Separator();
                for (int i = 2; i >= 0; --i) {
                    ImGui::MenuItem(plane_names[i], nullptr, &_crossSections[i].draw);
                }
                ImGui::EndPopup();
            }

            if (ImGui::BeginPopup("nematic_plane_selector")) {
                ImGui::TextDisabled("Active Nematic Overlays");
                ImGui::Separator();
                for (int i = 2; i >= 0; --i) {
                    bool& target = _widget.couplePlaneAndNematic ? _crossSections[i].draw : _crossSections[i].draw_nematic;
                    ImGui::MenuItem(plane_names[i], nullptr, &target);
                }
                ImGui::EndPopup();
            }

            if (_widget.couplePlaneAndNematic) {
                for (int i = 0; i < 3; ++i) {
                    _crossSections[i].draw_nematic = _crossSections[i].draw;
                }
            }

            ImGui::Separator();
            static bool s_flip_manipulator_rot = false;
            Matrix4 inversion_mat = Matrix4::rotationX(Rad(s_flip_manipulator_rot * M_PI));

            ImGui::AlignTextToFramePadding();
            ImGui::Text("Camera View:");
            ImGui::SameLine();
            if (ImGui::Button("XY")) {
                _manipulator->setTransformation(inversion_mat);
            }
            ImGui::SameLine();
            if (ImGui::Button("XZ")) {
                _manipulator->setTransformation(
                    inversion_mat * Matrix4::rotationZ(Rad(M_PI)) * Matrix4::rotationY(Rad(M_PI)) * Matrix4::rotationX(Rad(M_PI / 2.0f))
                );
            }
            ImGui::SameLine();
            if (ImGui::Button("YZ")) {
                _manipulator->setTransformation(
                    inversion_mat * Matrix4::rotationZ(Rad(-M_PI / 2.0f)) * Matrix4::rotationY(Rad(-M_PI / 2.0f))
                );
            }
            ImGui::SameLine();
            ImGui::Checkbox("Flip##Manipulator", &s_flip_manipulator_rot);
            ImGui::Spacing();
        }

        if (ImGui::CollapsingHeader("Analysis & Optical Toolboxes", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox("Enable POM", &_widget.POM);
            ImGui::SameLine();
            if (ImGui::Button("Refresh Visuals", ImVec2(-FLT_MIN, 0.f)) || _widget.updateImageFromLoad) {
                _widget.updateImage = true;
            }
            _widget.updateImageFromLoad = false;

            handleModificationWindow();
            handleInteractionEnergyWindow();
            handleNonlinearImagingWindow();
            handleZProfileWindow();
            handlePOMWindow();
            handlePreimageWindow();
            handleVortexKnotWindow();
            handleKnotInteractionWindow();
            handleMultiplaneWindow();
            handleLehmannWindow();
            ImGui::Spacing();
        }

        if (ImGui::CollapsingHeader("Solver & Relaxation", ImGuiTreeNodeFlags_DefaultOpen)) {
            auto relaxMap = Dataset::RelaxMap();
            dropDownMenu<Dataset::RelaxKind>("Method", data->relaxKind, relaxMap);

            if (ImGui::BeginTable("SolverParamsTable", 2, ImGuiTableFlags_SizingStretchSame)) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputInt("Cycle", &_widget.cycle);

                ImGui::TableNextColumn();
                Float relaxRate = data->rate;
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::InputFloat("Rate", &relaxRate, 0.0f, 0.0f, "%.4f")) {
                    data->rate = relaxRate;
                }
                ImGui::EndTable();
            }

            ImGui::RadioButton("Total Energy", &_widget.radioEn, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Energy Derivative", &_widget.radioEn, 1);

            _widget.relax = ImGui::Button("Relax Step", ImVec2(110.f, 0.f));
            ImGui::SameLine();
            ImGui::Checkbox("Auto-Relax", &_widget.autoRelax);
            ImGui::SameLine();
            ImGui::Checkbox("Hopf Charge", &_widget.sampleCharge);

            if (_widget.relax && _widget.autoRelax) {
                _widget.continuousRelax = true;
            }
            else if (!_widget.autoRelax) {
                _widget.continuousRelax = false;
            }

            if (_relaxFuture.second) {
                LC::scalar difference = 2.0 * _widget.energyErrorThreshold;
                if (_widget.energy_series.size() > 1) {
                    auto it_end = std::prev(_widget.energy_series.end());
                    auto it_prev = std::prev(it_end);
                    difference = *it_end - *it_prev;
                }

                ImGui::Spacing();
                ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.2f, 1.0f), "Status: Relaxing...");
                if (_widget.energy_series.size() > 1 && _widget.cycle > 0) {
                    ImGui::SameLine();
                    ImGui::TextDisabled("[dE/cycle = %.4e]", static_cast<float>(difference) / static_cast<float>(_widget.cycle));
                }
            }
        }

        ImGui::Separator();
        ImGui::TextDisabled("Frame: %.2f ms (%.1f FPS)",
            1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));

        ImGui::End();
        ImGui::PopStyleVar(4);
    }

    handleLCINFOWindow();

    if (_widget.showDemoWindow) {
        ImGui::SetNextWindowPos(ImVec2(100.f, 20.f), ImGuiCond_FirstUseEver);
        ImGui::ShowDemoWindow();
    }

    _imgui.updateApplicationCursor(*this);

    bool moving = true;
    if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::ArcBall)) {
        moving = _arcballCamera->updateTransformation();
    }

    updateBoxes();
    drawBoxes();

    updateSpheres();
    if (!_widget.hideHeliknotonInsertIndicator) {
        drawSpheres();
    }

    polyRenderer();

    sortObjects(_transparentNormalDrawables);
    if (_widget.drawSurfaces) {
        _camera->draw(_transparentNormalDrawables);
    }

    sortObjects(_transparentDrawables);
    _camera->draw(_transparentDrawables);

    const Matrix4 viewProj = _camera->cameraMatrix() * _manipulator->transformationMatrix();
    const Matrix4 projMat = _camera->projectionMatrix();

    for (int id = 0; id < 3; ++id) {
        if (_crossSections[id].draw_nematic && _crossSections[id].nematic) {
            auto& nem = _crossSections[id].nematic;
            nem->specular = _widget.specular;
            nem->diffuse = _widget.diffuse;
            nem->ambient = _widget.ambient;
            nem->selected_drawType = _widget.nematicDrawType;
            nem->Draw(viewProj, projMat);
        }
    }

    _multiplane_manager.Draw(viewProj, projMat, _widget);

    for (auto& knot : _processedVortexKnot) {
        if (knot.draw) {
            knot.points.Draw(viewProj, projMat);
        }
    }

    if (_widget.drawProfile && _zprofile) {
        _zprofile->Draw(viewProj, projMat);
    }

    guiRenderer();
    _imgui.drawFrame();

    swapBuffers();

    auto launchAsyncRelax = [&]() {
        _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, true);
        _relaxFuture.second = true;
        _widget.updateImage = true;
        };

    if (_widget.relax && !_widget.continuousRelax) {
        launchAsyncRelax();
    }

    if ((_widget.updateImage || _relaxFuture.second) && checkRelax()) {
        if (_widget.updateImage) {
            updateColor();
            computeEnergy();
            if (_widget.sampleCharge) {
                _widget.computed_topological_charge = LC::Math::ComputeHopfCharge(data->directors.get(), data->voxels);
            }
            _widget.updateImage = false;
        }
    }

    if (_widget.continuousRelax && checkRelax()) {
        launchAsyncRelax();
    }

    if (moving || _ioUpdate) {
        redraw();
    }
}

void Sandbox::updateBoxes(bool first) {
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    if (!first) arrayResize(_boxInstanceData, 0);
    
    Matrix4 camTransformationMatrix = _camera->cameraMatrix() * _manipulator->transformationMatrix();

    // Can just add setting in main window to hide...
    if (_widget.showCellBd)
        arrayAppend(_boxInstanceData, InPlaceInit,
            camTransformationMatrix *
            Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0], (float)data->cell_dims[1], (float)data->cell_dims[2] }), 0x00ffff_rgbf);

    if (_widget.showZProfileWindow && _zprofile) {

        arrayAppend(_boxInstanceData, InPlaceInit,
            camTransformationMatrix *
            Matrix4::translation(_zprofile->box.translation)
            * Matrix4::scaling(0.5f * _zprofile->box.dims),
            0x00ffff_rgbf
        );

    }

    if (_widget.showModificationWindow) {

        float dx = float(_widget.shrink_interval_end[0] - _widget.shrink_interval_begin[0]) / (data->voxels[0] - 1);
        float dy = float(_widget.shrink_interval_end[1] - _widget.shrink_interval_begin[1]) / (data->voxels[1] - 1);
        float dz = float(_widget.shrink_interval_end[2] - _widget.shrink_interval_begin[2]) / (data->voxels[2] - 1);

        float xb = 0.5f * (_widget.shrink_interval_begin[0] + _widget.shrink_interval_end[0]-2) / (data->voxels[0] - 1);
        float yb = 0.5f * (_widget.shrink_interval_begin[1] + _widget.shrink_interval_end[1]-2) / (data->voxels[1] - 1);
        float zb = 0.5f * (_widget.shrink_interval_begin[2] + _widget.shrink_interval_end[2]-2) / (data->voxels[2] - 1);

        arrayAppend(_boxInstanceData, InPlaceInit,
            camTransformationMatrix *
            Matrix4::translation(Vector3{ (- 0.5f + xb) * (float)data->cell_dims[0],
                (-0.5f + yb) * (float)data->cell_dims[1],
                (-0.5f + zb) * (float)data->cell_dims[2] }) *
            Matrix4::scaling(0.5f * Vector3{ (float)data->cell_dims[0] * dx,
                (float)data->cell_dims[1] * dy,
                (float)data->cell_dims[2] * dz }), 0x00ffff_rgbf);
    }

    
}

void Sandbox::updateSpheres(bool first) {
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    if (!first) arrayResize(_sphereInstanceData, 0);

    Matrix4 camTransformationMatrix = _camera->cameraMatrix() * _manipulator->transformationMatrix();

    if (_widget.showModificationWindow) {

        float xb = 1.f * (_widget.heliknotonInsertLocation[0] - 1) / (data->voxels[0] - 1);
        float yb = 1.f * (_widget.heliknotonInsertLocation[1] - 1) / (data->voxels[1] - 1);
        float zb = 1.f * (_widget.heliknotonInsertLocation[2] - 1) / (data->voxels[2] - 1);

        arrayAppend(_sphereInstanceData, InPlaceInit,
            camTransformationMatrix *
            Matrix4::translation(Vector3{ (-0.5f + xb) * (float)data->cell_dims[0],
                (-0.5f + yb) * (float)data->cell_dims[1],
                (-0.5f + zb) * (float)data->cell_dims[2] }) *
            Matrix4::scaling(0.1f * Vector3{ 1.f, 1.f, 1.f}), 0x00b800_rgbf);
    }
}

void Sandbox::drawBoxes() {
    _boxInstanceBuffer.setData(_boxInstanceData, GL::BufferUsage::DynamicDraw);
    _boxMesh.setInstanceCount(_boxInstanceData.size());
    _boxShader.setTransformationProjectionMatrix(_camera->projectionMatrix())
        .draw(_boxMesh);
}

void Sandbox::drawSpheres() {
    _sphereInstanceBuffer.setData(_sphereInstanceData, GL::BufferUsage::DynamicDraw);
    _sphereMesh.setInstanceCount(_sphereInstanceData.size());
    _boxShader.setTransformationProjectionMatrix(_camera->projectionMatrix())
        .draw(_sphereMesh);
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

    if ((event.key() == KeyEvent::Key::Up) || (event.key() == KeyEvent::Key::Down)) { // Zoom
        bool up_key = (event.key() == KeyEvent::Key::Up) ? 1 : 0;
        Float delta = .15f;

        // Fine focus
        if (event.modifiers() & KeyEvent::Modifier::Ctrl)
            delta *= .1f;

        if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::ArcBall)) {
            if (!up_key) delta *= -1.f;

            if (Magnum::Math::abs(delta) >= 1.0e-2f) _arcballCamera->zoom(delta);
        }
        if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::Group)) {

            Float cdelta = 1. - delta;

            /* Distance to origin */
            const Float distance = _cameraObject.transformation().translation().z();

            /* Move 15% of the distance back or forward */
            _cameraObject.translate(Vector3::zAxis(
                distance * (1.0f - (up_key ? 1 / cdelta : cdelta))));
        }
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

        // Call POM
        if (_widget.POM && ax == Axis::z) {
            POM();
            // Set nematic draw to false so that POM can be seen
            _crossSections[id].draw_nematic = false;
        }


        int xx = (ax == Axis::x) ? 1 : 0;
        int yy = (ax == Axis::y) ? 2 : xx + 1;

        std::size_t permutedSlice = data->voxels[xx];
        std::size_t permutedSlice2 = data->voxels[yy];

        int reduction_factor = _crossSections[id].nematic_reduction_factor;

        auto cross_idx = [permutedSlice](int i, int j) {
            return permutedSlice * j + i;
        };

        auto cross_idx_reduced = [permutedSlice, permutedSlice2, reduction_factor](int i, int j) {
            LC::scalar pos1 = (i / LC::scalar(permutedSlice - 1)); // [0, 1]
            LC::scalar pos2 = (j / LC::scalar(permutedSlice2 - 1)); // [0, 1]

            // Convert pos1 and pos2 to reduced space
            int reduced_slice1 = (permutedSlice + reduction_factor - 1) / reduction_factor;
            int reduced_slice2 = (permutedSlice2 + reduction_factor - 1) / reduction_factor;

            int i_reduced = ceil(pos1 * (reduced_slice1 - 1.1));
            int j_reduced = ceil(pos2 * (reduced_slice2 - 1.1));

            return reduced_slice1 * j_reduced + i_reduced;
        };

        int hvox;
        if (_widget.midplane) hvox = data->voxels[id] / 2;
        else hvox = _widget.iPlane[id];



        Vector3 scale((std::min)({ float(data->cell_dims[0]) / (data->voxels[0] - 1),
            float(data->cell_dims[1]) / (data->voxels[1] - 1),
            float(data->cell_dims[2]) / (data->voxels[2] - 1) }));

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
                phi = atan2(ny, nx);
                if (phi < 0) phi += 2.f * M_PI;

                std::size_t cidx = cross_idx(i, j);
                std::size_t cidx_red = cross_idx_reduced(i, j);
                Color3 director_color;

                if (!_widget.nonlinear) {
                    director_color = LC::Imaging::Colors::RungeSphere(theta, phi);
                }
                // Green 3 photon nonlinear imaging
                else {
                    if (_widget.nonlinCircular)
                        director_color = Color3::fromHsv({ Deg(120.0f), 1.0f, 1.0f - powf(nz, 6.0f) });
                    else {
                        director_color = Color3::fromHsv({ Deg(120.0f), 1.0f, powf(nx * cos(M_PI / 180. * _widget.nonlinTheta) + ny * sin(M_PI / 180. * _widget.nonlinTheta), 6.0f) });
                    }

                }
                if (_widget.S2colors == 1) _crossSections[id].nematic->polyInstanceData[cidx_red].color = director_color;
                else if (_widget.S2colors == 2) {
                    director_color = LC::Imaging::Colors::RungeSphere(M_PI / 2.f - abs(theta - M_PI/2.f), fmod(2.f * phi, 2.f * M_PI));
                    _crossSections[id].nematic->polyInstanceData[cidx_red].color = director_color;
                }
                else _crossSections[id].nematic->polyInstanceData[cidx_red].color = Color3{ 0.5f, 0.5f, 0.5f };
                if (!_widget.POM || ax != Axis::z) {
                    _crossSections[id].section.second->vertices[cidx].color = { director_color, alpha };
                }

                Vector3 translation{ 0.0f, 0.0f, 0.0f };
                translation[id] = data->cell_dims[id] * ((float)hvox / float(data->voxels[id] - 1) - 0.5f);
                translation[xx] = data->cell_dims[xx] / (data->voxels[xx] - 1) * i - data->cell_dims[xx] * .5f;
                translation[yy] = data->cell_dims[yy] / (data->voxels[yy] - 1) * j - data->cell_dims[yy] * .5f;

                _crossSections[id].section.second->vertices[cidx].position[id] = translation[id];
                float hPi = M_PI / 2.f;
                Matrix4 transformation = Matrix4::translation(translation)
                    * Matrix4::rotationZ(Rad{ -hPi + float(phi) })
                    * Matrix4::rotationX(Rad{ hPi - float(theta) })
                    * Matrix4::scaling(0.2f * scale * reduction_factor);

                // Map cidx to reduced space
                _crossSections[id].nematic->polyInstanceData[cidx_red].transformationMatrix = transformation;
            }
        }

        // Update sheet
        Trade::MeshData meshData = _crossSections[id].section.second->Data();
        _crossSections[id].section.second->vertexBuffer.setData(MeshTools::interleave(meshData.positions3DAsArray(), meshData.colorsAsArray()), GL::BufferUsage::DynamicDraw);
    }
}

void Sandbox::POM(std::string override_name) {
    if (!_solver) return;

    auto* data = static_cast<Dataset*>(_solver->GetDataPtr());
    if (!data || !data->directors) return;

    if (data->voxels[0] == 0 || data->voxels[1] == 0 || data->voxels[2] <= 1) {
        return;
    }

    int zSectionIdx = -1;
    for (int id = 0; id < 3; ++id) {
        if (_crossSections[id].axis == Axis::z) {
            zSectionIdx = id;
            break;
        }
    }

    if (zSectionIdx == -1 || !_crossSections[zSectionIdx].section.second) {
        return;
    }

    auto& targetSheet = *_crossSections[zSectionIdx].section.second;
    const std::size_t planeSize = static_cast<std::size_t>(data->voxels[0]) * data->voxels[1];

    if (targetSheet.vertices.size() < planeSize) {
        return;
    }

    const Float alpha = std::clamp(_widget.alpha, Float(0), Float(1));
    const LC::SIscalar pitch = _widget.pitch;
    const LC::scalar dop = data->cell_dims[2];

    if (data->lc_type != LC::FrankOseen::LC_TYPE::CUSTOM) {
        using OptConst = LC::FrankOseen::OpticalConstants;
        _pomImager.n0 = OptConst::LC(data->lc_type, OptConst::Constant::n_o).first;
        _pomImager.ne = OptConst::LC(data->lc_type, OptConst::Constant::n_e).first;
    }
    else {
        _pomImager.n0 = data->n0;
        _pomImager.ne = data->ne;
    }

    _pomImager.dop = dop;
    _pomImager.thickness = pitch.first * (dop + _pomImager.additional_layers);
    _pomImager.dz = (pitch.first * dop * 1e-6) / static_cast<double>(data->voxels[2] - 1);

    auto assignVertexColor = [](void* rawVertices, const std::array<float, 4>& color, std::size_t idx) {
        auto* vertices = static_cast<Magnum::Containers::Array<LC::DynamicColorSheet::Vertex>*>(rawVertices);
        for (int channel = 0; channel < 4; ++channel) {
            (*vertices)[idx].color[channel] = color[channel];
        }
        };

    _pomImager.Compute(
        data->directors.get(),
        data->voxels,
        static_cast<void*>(&targetSheet.vertices),
        assignVertexColor,
        alpha
    );

    if (_widget.savePOM) {
        using namespace LC::Imaging;
        auto plane = std::make_unique<ImageSeries::COLOR[]>(planeSize);

        auto toByte = [](float val) -> unsigned char {
            return static_cast<unsigned char>(std::round(std::clamp(val, 0.0f, 1.0f) * 255.0f));
            };

        const unsigned char alphaByte = toByte(static_cast<float>(alpha));

        for (std::size_t idx = 0; idx < planeSize; ++idx) {
            const auto& vertexColor = targetSheet.vertices[idx].color;
            plane[idx].R = toByte(vertexColor[0]);
            plane[idx].G = toByte(vertexColor[1]);
            plane[idx].B = toByte(vertexColor[2]);
            plane[idx].A = alphaByte;
        }

        const std::string rawSource = override_name.empty() ? getLoadFile() : override_name;
        std::string baseName = std::filesystem::path(rawSource).stem().string();

        if (baseName.empty()) {
            baseName = "default";
        }

        std::filesystem::path outputPath = std::filesystem::path(_widget.savePOM_loc) / baseName;
        _image_series.write_file = outputPath.generic_string();

        LC_INFO("Saving POM as <{0}>", _image_series.write_file.c_str());
        _image_series.GenerateAndWriteFrame(plane.get(), data->voxels[0], data->voxels[1]);
    }

    Trade::MeshData meshData = targetSheet.Data();
    targetSheet.vertexBuffer.setData(
        MeshTools::interleave(meshData.positions3DAsArray(), meshData.colorsAsArray()),
        GL::BufferUsage::DynamicDraw
    );
}

void Sandbox::initVisuals() {
    if (!_solver) return;

    auto* data = static_cast<Dataset*>(_solver->GetDataPtr());
    if (!data) return;

    const std::array<int, 3> voxels = data->voxels;
    const std::array<LC::scalar, 3> cdims = data->cell_dims;

    if (voxels[0] <= 1 || voxels[1] <= 1 || voxels[2] <= 1) {
        return;
    }

    _widget.shrink_interval_end = voxels;
    _widget.shrink_interval_begin = { 1, 1, 1 };
    _widget.radioZProfileIndex = -1;
    _widget.heliknotonInsertLocation = { voxels[0] / 2, voxels[1] / 2, voxels[2] / 2 };

    _zprofile = {};
    _quaternion_plane = {};
    if (_pionPreimage) {
        _pionPreimage->draw = false;
    }
    _pionPreimage = {};
    _baryon_density.reset();
    _quaternion_field.reset();

    _vortex_line.components.clear();
    _vortexKnot.clear();
    _processedVortexKnot.clear();

    for (int d = 0; d < 3; ++d) {
        _widget.celldims[d] = cdims[d];
        _widget.boundaries[d] = data->bc[d];
    }

    _preimageManipulator.reset();
    _vortexManipulator.reset();
    _processedVortexManipulator.reset();

    _manipulator = std::make_unique<LC::Drawable::Object3D>(&_scene);
    _preimageManipulator = std::make_unique<LC::Drawable::Object3D>(_manipulator.get());
    _vortexManipulator = std::make_unique<LC::Drawable::Object3D>(_manipulator.get());
    _processedVortexManipulator = std::make_unique<LC::Drawable::Object3D>(_manipulator.get());

    _cameraObject.setTransformation(
        Matrix4::translation(Vector3::zAxis(2.2f * (std::max)(cdims[0], cdims[1])))
    );
    _manipulator->setTransformation(Matrix4{});

    _crossSections = Containers::Array<CrossX>{ 3 };
    _multiplane_manager.CleanAll();

    const std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* rawPos, std::size_t idx) {
        return static_cast<Magnum::Vector3*>(rawPos)[idx];
        };

    for (int id = 0; id < 3; ++id) {
        const Axis ax = static_cast<Axis>(id);
        const int i = (ax == Axis::x) ? 1 : 0;
        const int j = (ax == Axis::y) ? 2 : i + 1;

        if (ax != Axis::z) {
            _crossSections[id].draw = false;
            _crossSections[id].draw_nematic = false;
        }

        const int next = (std::max)(1, _crossSections[id].nematic_reduction_factor);
        _crossSections[id].nematic = LC::NematicArray{};

        const int vox_red_x = (voxels[i] + next - 1) / next;
        const int vox_red_y = (voxels[j] + next - 1) / next;
        const std::size_t size_reduced = static_cast<std::size_t>(vox_red_x) * vox_red_y;

        auto positions = std::make_unique<Vector3[]>(size_reduced);
        const float dx = static_cast<float>(cdims[i]) / static_cast<float>(voxels[i] - 1);
        const float dy = static_cast<float>(cdims[j]) / static_cast<float>(voxels[j] - 1);
        const float half_x = 0.5f * static_cast<float>(cdims[i]);
        const float half_y = 0.5f * static_cast<float>(cdims[j]);

        for (int x = 0, x_red = 0; x < voxels[i]; x += next, ++x_red) {
            for (int y = 0, y_red = 0; y < voxels[j]; y += next, ++y_red) {
                const std::size_t flat_idx = static_cast<std::size_t>(x_red) + static_cast<std::size_t>(y_red) * vox_red_x;
                Vector3 pos{ 0.0f, 0.0f, 0.0f };
                pos[i] = -half_x + x * dx;
                pos[j] = -half_y + y * dy;
                positions[flat_idx] = pos;
            }
        }

        _crossSections[id].nematic->scale *= next;
        _crossSections[id].nematic->Init(positions.get(), access, size_reduced);
        _crossSections[id].axis = ax;

        _crossSections[id].section.second = std::make_unique<LC::DynamicColorSheet>();
        _crossSections[id].section.second->Set(voxels[i], voxels[j], cdims[i], cdims[j]);
        _crossSections[id].section.second->Init(CrossX::Position[id], 0.0f);
        _crossSections[id].section.first = _crossSections[id].section.second->Mesh();

        new LC::Drawable::TransparentFlatDrawable{
            *_manipulator,
            _transparentShader,
            *_crossSections[id].section.first,
            _crossSections[id].draw,
            Vector3{ 0.0f, 0.0f, 0.0f },
            _transparentDrawables
        };
    }
}

void Sandbox::imageMenu() {
    if (ImGui::BeginMenuBar()) {
        if (ImGui::BeginMenu("Tools")) {
            if (ImGui::MenuItem("Isosurfaces")) { _widget.showPreimageSettings = true; }
            if (ImGui::MenuItem("Vortex Knot")) { _widget.showVortexKnotSettings = true; }
            if (ImGui::MenuItem("Nonlinear Imaging")) { _widget.showNonlinearSettings = true; }
            if (ImGui::MenuItem("Z-Profile")) { _widget.showZProfileWindow = true; }
            if (ImGui::MenuItem("Knot Interaction")) { _widget.knot_interaction_handle.showWindow = true; }
            if (ImGui::MenuItem("Interaction Energy")) { _widget.showInteractionEnergyWindow = true; }
            if (ImGui::MenuItem("Multi-plane")) { _widget.multiplane_window = true; }
            if (ImGui::MenuItem("Draw Lehmann Clusters")) { _widget.lehmann_tool_window = true; }

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Settings")) {
            if (ImGui::MenuItem("Visual properties")) { _widget.showVisualProperties = true; }
            if (ImGui::MenuItem("POM")) { _widget.showPOMSettings = true; }

            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }
}

// Helper to preview a wavelength (nm) as an approximate display RGB color in ImGui
static ImVec4 wavelengthToDisplayRGB(float lambda_nm) {
    float r = 0.0f, g = 0.0f, b = 0.0f;
    if (lambda_nm >= 380.0f && lambda_nm < 440.0f) {
        r = -(lambda_nm - 440.0f) / (440.0f - 380.0f);
        b = 1.0f;
    }
    else if (lambda_nm >= 440.0f && lambda_nm < 490.0f) {
        g = (lambda_nm - 440.0f) / (490.0f - 440.0f);
        b = 1.0f;
    }
    else if (lambda_nm >= 490.0f && lambda_nm < 510.0f) {
        g = 1.0f;
        b = -(lambda_nm - 510.0f) / (510.0f - 490.0f);
    }
    else if (lambda_nm >= 510.0f && lambda_nm < 580.0f) {
        r = (lambda_nm - 510.0f) / (580.0f - 510.0f);
        g = 1.0f;
    }
    else if (lambda_nm >= 580.0f && lambda_nm < 645.0f) {
        r = 1.0f;
        g = -(lambda_nm - 645.0f) / (645.0f - 580.0f);
    }
    else if (lambda_nm >= 645.0f && lambda_nm <= 780.0f) {
        r = 1.0f;
    }
    return ImVec4(r, g, b, 1.0f);
}

void Sandbox::handlePOMWindow() {
    if (!_widget.showPOMSettings) return;

    ImGui::SetNextWindowSize(ImVec2(420, 540), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Configure POM", &_widget.showPOMSettings)) {
        ImGui::End();
        return;
    }

    if (ImGui::CollapsingHeader("Sample & Geometry", ImGuiTreeNodeFlags_DefaultOpen)) {
        float pitch = static_cast<float>(_widget.pitch.first);
        if (ImGui::DragFloat("Pitch (um)", &pitch, 0.05f, 0.01f, 1000.0f, "%.3f")) {
            _widget.pitch.first = (std::max)(0.001f, pitch);
        }

        ImGui::SliderFloat("z-depth ratio", &_pomImager.z_scan_ratio, 0.0f, 1.0f, "%.3f");

        if (ImGui::InputInt("+ layers (p/2)", &_pomImager.additional_layers)) {
            _pomImager.additional_layers = (std::max)(0, _pomImager.additional_layers);
        }
    }

    if (ImGui::CollapsingHeader("Optical Train", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Polarizer", &_pomImager.polarizer);
        ImGui::SameLine(140.0f);
        ImGui::Checkbox("Analyzer", &_pomImager.analyzer);

        ImGui::SliderFloat("Polarizer Angle", &_pomImager.polarizerAngle, -180.0f, 180.0f, "%.1f deg");
        ImGui::DragFloat("Gamma", &_pomImager.gamma, 0.02f, 0.1f, 5.0f, "%.2f");

        using Waveplate = LC::Imaging::UniformGrid::POM::Waveplate;
        std::map<Waveplate, std::string> waveplateMap{
            { Waveplate::None, "None" },
            { Waveplate::Full530nm, "530 nm Full-Wave" }
        };
        dropDownMenu<Waveplate>("Waveplate", _pomImager.waveplate, waveplateMap);
    }

    if (ImGui::CollapsingHeader("Light Spectrum", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Button("+ Add Wavelength")) {
            _pomImager.lights.push_back({ 550.0f, 1.0f });
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset RGB (3)")) {
            _pomImager.lights = { { 650.0f, 1.0f },{ 550.0f, 0.6f },{ 450.0f, 0.2f } };
        }
        ImGui::SameLine();
        if (ImGui::Button("Broadband (9)")) {
            _pomImager.lights = {
                { 420.0f, 0.25f },{ 450.0f, 0.60f },{ 480.0f, 0.80f },
                { 510.0f, 0.90f },{ 540.0f, 1.00f },{ 570.0f, 0.95f },
                { 600.0f, 0.85f },{ 630.0f, 0.70f },{ 660.0f, 0.45f }
            };
        }

        const ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersInnerV |
            ImGuiTableFlags_RowBg |
            ImGuiTableFlags_SizingStretchProp;

        if (ImGui::BeginTable("LightsTable", 4, tableFlags)) {
            ImGui::TableSetupColumn("##Swatch", ImGuiTableColumnFlags_WidthFixed, 22.0f);
            ImGui::TableSetupColumn("Wavelength (nm)", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("Rel. Intensity", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("##Delete", ImGuiTableColumnFlags_WidthFixed, 28.0f);
            ImGui::TableHeadersRow();

            int removeIdx = -1;
            for (int i = 0; i < static_cast<int>(_pomImager.lights.size()); ++i) {
                ImGui::PushID(i);
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                ImVec4 swatch = wavelengthToDisplayRGB(_pomImager.lights[i].wavelength);
                ImGui::ColorButton("##wl_color", swatch, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoDragDrop, ImVec2(18, 18));

                ImGui::TableSetColumnIndex(1);
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragFloat("##wl", &_pomImager.lights[i].wavelength, 1.0f, 380.0f, 780.0f, "%.1f nm");

                ImGui::TableSetColumnIndex(2);
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragFloat("##int", &_pomImager.lights[i].intensity, 0.01f, 0.0f, 10.0f, "%.2f");

                ImGui::TableSetColumnIndex(3);
                if (_pomImager.lights.size() > 1) {
                    if (ImGui::SmallButton("X")) {
                        removeIdx = i;
                    }
                }

                ImGui::PopID();
            }

            if (removeIdx >= 0) {
                _pomImager.lights.erase(_pomImager.lights.begin() + removeIdx);
            }

            ImGui::EndTable();
        }
    }

    if (ImGui::CollapsingHeader("Export", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Save POM on Compute", &_widget.savePOM);
        if (_widget.savePOM) {
            if (ImGui::Button("Choose Directory...")) {
                auto res = pfd::select_folder("POM Location", LCLAB2_ROOT_PATH, pfd::opt::force_path);
                while (!res.ready(1000)) {
                    LC_INFO("Waiting for user...");
                }
                if (res.ready() && !res.result().empty()) {
                    _widget.savePOM_loc = res.result();
                }
            }
            ImGui::SameLine();
            ImGui::TextWrapped("%s", _widget.savePOM_loc.c_str());
        }
    }
    ImGui::End();
}

void Sandbox::handleVisualPropertiesWindow() {
    if (!_widget.showVisualProperties)
        return;

    ImGui::Begin("Visual properties", &_widget.showVisualProperties);

    // Background color
    ImGui::ColorEdit4("Background Color##BackgroundColor", &_widget.clearColor[0], ImGuiColorEditFlags_NoInputs);
    ImGui::Separator();
    ImGui::Text("Director color scheme");
    ImGui::Separator();
    ImGui::RadioButton("Gray", &_widget.S2colors, 0);
    ImGui::SameLine();
    ImGui::RadioButton("S2", &_widget.S2colors, 1);
    ImGui::SameLine();
    ImGui::RadioButton("S2/Z2", &_widget.S2colors, 2);

    ImGui::Separator();
    ImGui::PushItemWidth(75.f);
    ImGui::InputFloat("Specular", &_widget.specular);
    ImGui::SameLine();
    ImGui::InputFloat("Diffuse", &_widget.diffuse);
    ImGui::SameLine();
    ImGui::InputFloat("Ambient", &_widget.ambient);
    ImGui::PopItemWidth();

    ImGui::SliderFloat("Plane alpha", &_widget.alpha, 0.0f, 1.0f);
    ImGui::Checkbox("Show cell boundary", &_widget.showCellBd);

    dropDownMenu<LC::NematicArray::DrawType>("Nematic draw type", _widget.nematicDrawType, LC::NematicArray::DrawMap());

    ImGui::End();
}

void Sandbox::handlePreimageWindow() {
    if (_widget.showPreimageSettings) {
        Dataset* data = (Dataset*)(_solver->GetDataPtr());

        ImGui::Begin("Configure Preimages", &_widget.showPreimageSettings);

        if (ImGui::CollapsingHeader("Global Properties")) {
            ImGui::PushItemWidth(100.0f);

            ImGui::SliderFloat("Preimage alpha", &_widget.preimage_alpha, 0.0f, 1.0f);
            ImGui::SameLine();
            ImGui::Checkbox("Global preimage alpha", &_widget.global_preimage_alpha);
            ImGui::InputInt("Smoothing iterations", &_widget.smoothingIterations);
            ImGui::SliderFloat("Laplacian Smoothing coeff", &_widget.smoothingAlpha, 0.25f, 200.f);
            ImGui::SliderFloat("Taubin Smoothing contraction coeff", &_widget.smoothingLambda, 0.0f, 0.5f);
            ImGui::SliderFloat("Taubin Smoothing expansion coeff", &_widget.smoothingMu, -0.5f, 0.0f);
            ImGui::TextColored({0.f, 1.f, 0.f, 1.f}, "Smoothing type");
            ImGui::RadioButton("Explicit", &_widget.smoothingType, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Iterative", &_widget.smoothingType, 1);
            ImGui::SameLine();
            ImGui::RadioButton("Implicit", &_widget.smoothingType, 2);

            ImGui::PopItemWidth();

            int maxVox = (std::max)({ data->voxels[0], data->voxels[1], data->voxels[2] });
            ImGui::PushItemWidth(250.f);
            ImGui::SliderInt3("Start Indices", &_widget.shrink_interval_begin[0], 0, maxVox);
            ImGui::SliderInt3("End Indices", &_widget.shrink_interval_end[0], 0, maxVox);
            ImGui::PopItemWidth();

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

            if (!_widget.showModificationWindow) {
                updateBoxes();
                drawBoxes();
            }
            ImGui::Separator();
        }

        ImGui::TextColored({ 1.f, 1.f, 0.f, 1.f }, "Nematic preimages");

        ImGui::PushItemWidth(100.0f);

        ImGui::SliderInt("theta", &_widget.ptheta, 0, 180);
        ImGui::SameLine();
        ImGui::SliderInt("phi", &_widget.pphi, 0, 360);
        ImGui::SameLine();
        ImGui::SliderFloat("isovalue##widget", &_widget.isoLevel, 0.0f, 0.3f);

        ImGui::PopItemWidth();

        ImGui::SameLine();

        if (ImGui::Button("Add preimage")) {
            // Append a preimage to the list

            // Make sure the preimage doesn't already exist!
            bool found = false;
            for (auto& p : _nematicPreimages) {

                if (p.theta == _widget.ptheta && p.phi == _widget.pphi) {
                    found = true;
                    p.isoLevel = _widget.isoLevel;

                }
            }
            if (!found)
                _nematicPreimages.emplace_back((float)_widget.ptheta, (float)_widget.pphi, _widget.isoLevel);
        }

        // Show tab bar to manage preimages
        nematicPreimageManager();

        if (ImGui::Button("Create Vortex Shell")) {
            // Create a new vortex knot
            _vortexShell = VortexShell{};
            _vortexShell->isoLevel = _widget.isoLevel;
        }

        ImGui::SameLine();

        if (ImGui::Button("Remove Vortex Shell")) {
            if (_vortexShell)
                _vortexShell->draw = false;
            _vortexShell = {};
        }

        if (_vortexShell) {
            ImGui::SameLine();
            ImGui::Checkbox("No cull", &_vortexShell->noCull);
            ImGui::SameLine();
            ImGui::Checkbox("Draw shell", &_vortexShell->draw);
            ImGui::SameLine();
            ImGui::Checkbox("Domain inversion", &_vortexShell->invertNormals);
            ImGui::PushItemWidth(100.0f);
            ImGui::SliderFloat("Shell alpha", &_vortexShell->alpha, 0.0f, 1.0f);
            ImGui::SameLine();
            ImGui::SliderFloat("Shell isovalue", &_vortexShell->isoLevel, 0.0f, 1.0f);
            ImGui::SliderFloat("Shell query value", & _vortexShell->queryValue, -1.f, 1.f);
            ImGui::RadioButton("Tilt (theta)", &_widget.chiColorScheme, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Tilt (chi_x)", &_widget.chiColorScheme, 1);
            ImGui::SameLine();
            ImGui::RadioButton("Quaternion", &_widget.chiColorScheme, 2);
            ImGui::PopItemWidth();
        }

        static bool generated_preimages = false;
        bool gen_surfaces = ImGui::Button("Generate surfaces");
        
        if (gen_surfaces || _widget.generateKnots) {
            if (gen_surfaces) {
                generated_preimages = true;
            }
            generateIsosurface();
        }

        if (_nematicPreimages.size() && generated_preimages) {
            // Show settings parametrizing conversion from surface to graph

            const size_t mod_pimage_graph_file_loc_sz = 128;
            static int graph_min_comp_size = 5;
            static int graph_knn = 4;
            static float exclusion_scale = 1.5f;
            static char pimage_fdir[mod_pimage_graph_file_loc_sz] = "../data/";
            static char fname_pimage_graph_export[mod_pimage_graph_file_loc_sz] = "pknot";

            ImGui::Text("Preimage to graph parameters");
            ImGui::Separator();
            ImGui::InputText("Base file directory##pimage", &pimage_fdir[0], mod_pimage_graph_file_loc_sz);
            ImGui::InputText("Base file name##pimage",&fname_pimage_graph_export[0], mod_pimage_graph_file_loc_sz);
            ImGui::InputInt("Minimum component size##pimage", &graph_min_comp_size);
            ImGui::InputFloat("Exclusion scale (s: 2*s*dr)##pimage", &exclusion_scale);
            ImGui::InputInt("knn##pimage", &graph_knn);

            if (ImGui::Button("Convert preimage(s) to graph(s)##pimage")) {
                // CODE: Nematic preimages is a good vector: convert to knot form and export

                // Define the rbf (for k nearest neighbors)
                LC::Math::poly_spline<LC::scalar> rbf;
                LC::Math::Metric<LC::scalar> metric;
                metric.SetBCS(LC::Math::PBC({ false, false, false }));
                metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[1]);
                auto vox = data->voxels;
                std::size_t slice = vox[0] * vox[1];
                std::size_t vol = slice * vox[2];

                // Translate each component into sets of points
                LC::scalar dx = data->cell_dims[0] / (vox[0] - 1);
                LC::scalar dy = data->cell_dims[1] / (vox[1] - 1);
                LC::scalar dz = data->cell_dims[2] / (vox[2] - 1);
                LC::scalar dr = sqrt(dx * dx + dy * dy + dz * dz);
                // Separation distance of two points
                LC::scalar min_point_dist = 2. * exclusion_scale * dr;


                for (auto& pimage : _nematicPreimages) {

                    std::vector<Eigen::Vector3d> distinct_points;
                    std::vector<int> npoints_per_component;
                    std::vector<Eigen::Vector3d> points;
                    
                    for (const auto& vertex : pimage.surface.vertices) {

                        Eigen::Vector3d potential_pos(vertex.position[0], vertex.position[1], vertex.position[2]);

                        // Check if position is too close to existing points
                        bool valid_pt = true;
                        for (const auto& pt : distinct_points) {
                            LC::scalar dist = (potential_pos - pt).norm();
                            if (dist < min_point_dist) {
                                valid_pt = false;
                                break;
                            }
                        }

                        // Add to the distinct points
                        if (valid_pt) {
                            distinct_points.emplace_back(potential_pos);
                        }
                    }

                    // See how many distinct points
                    LC_INFO("{0} distinct points detected", distinct_points.size());

                    LC::Math::Metric<LC::scalar> metric;
                    metric.Bcs = { false, false, false };
                    metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);
                    std::unique_ptr<LC::Math::rbf<LC::scalar>> RBF = std::unique_ptr<LC::Math::poly_spline<LC::scalar>>(new LC::Math::poly_spline<LC::scalar>);
                    std::unique_ptr<std::size_t[]> neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[distinct_points.size() * graph_knn]);
                    // Convert distinct points into proper data format for knn
                    std::unique_ptr<LC::scalar[]> positions(new LC::scalar[distinct_points.size() * 3]);
                    std::unique_ptr<std::size_t[]> query_domain(new std::size_t[distinct_points.size()]);


                    for (int i = 0; i < distinct_points.size(); i++) {
                        // Entire query domain for the knot
                        query_domain[i] = i;
                        for (int d = 0; d < 3; d++)
                            positions[i + d * distinct_points.size()] = (distinct_points[i])(d);
                    }

                    // Find nearest neighbors
                    LC::Algorithm::knn_c(positions.get(), distinct_points.size(),
                        query_domain.get(), distinct_points.size(),
                        metric, graph_knn, (LC::scalar*)0, neighbors.get());

                    // If a point is isolated, that is d(p,q) > 2*separation for all q != p, throw the point out
                    // Do this for all isolated points and recompute knn

                    std::vector<uint> discard_list;

                    for (auto i = 0; i < distinct_points.size(); i++) {
                        int mapoffset = distinct_points.size();
                        std::size_t n1 = neighbors[mapoffset + i];

                        auto pt = distinct_points[i];
                        auto ptn = distinct_points[n1];

                        LC::scalar dist = (pt - ptn).norm();

                        if (dist > 2. * min_point_dist) {
                            discard_list.emplace_back(i);
                        }
                    }

                    if (discard_list.size())
                        LC_INFO("Discarding {0} points", discard_list.size());

                    // Parse the list again
                    std::vector<Eigen::Vector3d> parsed_points;

                    for (int i = 0; i < distinct_points.size(); i++) {
                        bool remove_pt = false;
                        for (auto id : discard_list)
                            if (id == i)
                                remove_pt = true;

                        if (!remove_pt)
                            parsed_points.emplace_back(distinct_points[i]);
                    }

                    if (discard_list.size() > 0) {
                        // Refill position data
                        for (int i = 0; i < parsed_points.size(); i++) {
                            for (int d = 0; d < 3; d++)
                                positions[i + d * parsed_points.size()] = (parsed_points[i])(d);
                        }

                        // Recompute neighbors (No need to alter array sizes, there will just be surplus allocated space at the end)
                        LC::Algorithm::knn_c(positions.get(), parsed_points.size(),
                            query_domain.get(), parsed_points.size(),
                            metric, graph_knn, (LC::scalar*)0, neighbors.get());
                    }

                    // Create a graph with the remaining vertices
                    // Add edges
                    Graph graph(parsed_points.size());
                    int mapoffset = parsed_points.size();

                    for (int i = 0; i < parsed_points.size(); i++) {
                        // Compute the dist between each k nearest neighbors
                        // if the distance is smaller than the minimum required distance, add those neighbors as edges in the graph
                        for (int k = 1; k < graph_knn; k++) {
                            std::size_t nb_k = neighbors[k * mapoffset + i];

                            Eigen::Vector3d pt = parsed_points[i];
                            Eigen::Vector3d pt_nb_k = parsed_points[nb_k];

                            Eigen::Vector3d tk = pt - pt_nb_k;
                            LC::scalar diff = tk.norm();

                            if (diff < 2. * min_point_dist) {
                                graph.addEdge(i, nb_k);
                            }
                        }
                    }

                    if (parsed_points.size()) {
                        // Save the data
                        std::string spimage_fdir = std::string(pimage_fdir);
                        if (spimage_fdir.size()) {
                            // Make sure directory is valid
                            if (spimage_fdir[spimage_fdir.size() - 1] != '/' && spimage_fdir[spimage_fdir.size() - 1] != '\\') {
                                spimage_fdir += '/';
                            }
                        }

                        std::string pimage_classifier = "_theta_" + std::to_string((int)pimage.theta) + std::string("_phi_") + std::to_string((int)pimage.phi);
                        std::string pimage_fname = spimage_fdir + std::string(fname_pimage_graph_export) + pimage_classifier;

                        std::ofstream ofile_graph((pimage_fname + ".txt").c_str(), std::ios::out);

                        if (ofile_graph.is_open()) {                        
                            // 1. Save the points
                            // Save number of points
                            ofile_graph << parsed_points.size() << std::endl;
                            // Save the points
                            for (int i = 0; i < parsed_points.size(); i++) {
                                ofile_graph << (float)parsed_points[i].x() * 20.f << " "
                                    << (float)parsed_points[i].y() * 20.f << " "
                                    << (float)parsed_points[i].z() * 20.f << std::endl;
                            }

                            // 2. Save the adjacency list
                            auto adjList = graph.adjacencyList();

                            // Save the elements
                            for (int i = 0; i < parsed_points.size(); i++) {
                                std::list<uint>& aL_i = adjList[i];
                                std::size_t sz_i = aL_i.size();
                                std::size_t iter = 0;
                                for (auto elem : aL_i) {
                                    ofile_graph << elem;
                                    if (iter++ != sz_i - 1)
                                        ofile_graph << " "; // Add space if not the last element of list
                                }
                                if (i != parsed_points.size() - 1)
                                    ofile_graph << std::endl;
                            }
                        

                            // Close the file
                            ofile_graph.close();
                        }
                    }
                }
            }
        }
        else if (_nematicPreimages.size() == 0) {
            // Preimages got reset
            generated_preimages = false;
        }

        ImGui::Checkbox("Draw Surfaces", &_widget.drawSurfaces);

        ImGui::End();

    }
}

void Sandbox::handleKnotInteractionWindow() {
    if (_widget.knot_interaction_handle.GUI()) {
        Dataset* data = (Dataset*)(_solver->GetDataPtr());
        auto solver = (FOFDSolver*)_solver.get();

        if (_widget.knot_interaction_handle.Dispatch()) {

            // Check if file directory exists, if not make it!
            LC_MKDIR(_widget.knot_interaction_handle.PathToFile());

            // Begin the interaction routine
            std::vector<uint32_t> full_region;
            Eigen::Vector3d disp(0., 0., 0.);
            std::array<LC::scalar, 3> CELL = data->cell_dims;
            std::array<int, 3> VOX = data->voxels;

            LC::scalar dx = CELL[0] / (VOX[0] - 1);
            LC::scalar dy = CELL[1] / (VOX[1] - 1);
            LC::scalar dz = CELL[2] / (VOX[2] - 1);

            std::vector<Eigen::Vector3d> translations;

            // 1 Create the initial conditons (determined by settings)
            if (_widget.knot_interaction_handle.useInitialConditions) {
                // Create initial conditions

                /*
                    TEMPORARY INITIAL CONDITIONS
                */
                int NODE_DENSITY = _widget.knot_interaction_handle.node_density;
                CELL = _widget.knot_interaction_handle.GetCell();
                for (int i = 0; i < 3; i++) VOX[i] = NODE_DENSITY * CELL[i];

                data->directors = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * VOX[0] * VOX[1] * VOX[2]]);
                std::unique_ptr<LC::scalar[]> nn_data_temp(new LC::scalar[3 * VOX[0] * VOX[1] * VOX[2]]);
                data->voltage = std::unique_ptr<LC::scalar[]>(new LC::scalar[VOX[0] * VOX[1] * VOX[2]]);
                data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[VOX[0] * VOX[1] * VOX[2]]);

                FOFDSolver::Tensor4 nn(data->directors.get(), VOX[0], VOX[1], VOX[2], 3);
                FOFDSolver::Tensor3 vv(data->voltage.get(), VOX[0], VOX[1], VOX[2]);
                FOFDSolver::Tensor4 nn_temp(nn_data_temp.get(), VOX[0], VOX[1], VOX[2], 3);
                // Change voxels and cell dims
                data->voxels = VOX;
                data->cell_dims = CELL;

                // Modified voxels, make sure to update graphics

                dx = CELL[0] / (VOX[0] - 1);
                dy = CELL[1] / (VOX[1] - 1);
                dz = CELL[2] / (VOX[2] - 1);

                // Generate helical field
                auto helical = LC::Math::Planar(2, 1);

                auto clear_heliknotons = [&]() {
                    // Initialize to the uniform helical bg
                    for (int i = 0; i < VOX[0]; i++) {
                        for (int j = 0; j < VOX[1]; j++) {
                            for (int k = 0; k < VOX[2]; k++) {
                                LC::scalar z = -0.5 * CELL[2] + k * dz;
                                auto bg = helical(0., 0., z);

                                if (!_widget.knot_interaction_handle.switching_data.empty())
                                    vv(i, j, k) = z * _widget.knot_interaction_handle.switching_data[0].efield;
                                else
                                    vv(i, j, k) = 0.;

                                for (int d = 0; d < 3; d++) {
                                    nn(i, j, k, d) = bg[d];
                                }
                            }
                        }
                    }
                };


                auto embed_heliknoton = [&](Eigen::Vector3d translation_vector) {
                    LC::scalar phi0 = 2. * M_PI * data->chirality * translation_vector[2];
                    auto heli_vfield = LC::Math::Heliknoton(1, data->cell_dims, 1., 1., { 0.,0.,0. }, phi0, false);

                    for (int i = 0; i < VOX[0]; i++)
                        for (int j = 0; j < VOX[1]; j++)
                            for (int k = 0; k < VOX[2]; k++) {
                                // Current coords
                                Eigen::Vector3d pos;
                                pos[0] = data->cell_dims[0] * ((LC::scalar)i / (data->voxels[0] - 1) - 0.5);
                                pos[1] = data->cell_dims[1] * ((LC::scalar)j / (data->voxels[1] - 1) - 0.5);
                                pos[2] = data->cell_dims[2] * ((LC::scalar)k / (data->voxels[2] - 1) - 0.5);

                                pos = pos - translation_vector;

                                auto res = heli_vfield(pos[0], pos[1], pos[2]);
                                if (res[0] != 0. || res[1] != 0. || res[2] != 0.) {
                                    nn(i, j, k, 0) = res[0];
                                    nn(i, j, k, 1) = res[1];
                                    nn(i, j, k, 2) = res[2];
                                }
                            }
                };

                if (_widget.knot_interaction_handle.useInitialConditions == 1) {
                    // Define the heliknoton displacement vector using GUI specifications
                    LC::scalar r = 0.5 * _widget.knot_interaction_handle.seperation;
                    LC::scalar theta0 = M_PI / 180.f * _widget.knot_interaction_handle.theta0;
                    LC::scalar phi0 = M_PI / 180.f * _widget.knot_interaction_handle.phi0;


                    disp[0] = r * sin(theta0) * cos(phi0);
                    disp[1] = r * sin(theta0) * sin(phi0);
                    disp[2] = r * cos(theta0);


                    // Set the background field to uniform
                    clear_heliknotons();

                    translations.push_back(disp);
                    translations.push_back(-disp);

                    // Embed the heliknotons
                    embed_heliknoton(disp);
                    embed_heliknoton(-disp);
                }
                else if (_widget.knot_interaction_handle.useInitialConditions == 2) {
                    // Set the background field to uniform
                    clear_heliknotons();

                    for (auto& pos : _widget.knot_interaction_handle.pos_gui.position_array) {
                        translations.emplace_back(Eigen::Vector3d(pos[0], pos[1], pos[2]));
                        embed_heliknoton(Eigen::Vector3d(pos[0], pos[1], pos[2]));
                    }

                }
                else if (_widget.knot_interaction_handle.useInitialConditions == 3) {
                    // Define the heliknoton displacement vector using GUI specifications
                    LC::scalar r = 0.5 * _widget.knot_interaction_handle.seperation;
                    LC::scalar theta0 = M_PI / 180.f * _widget.knot_interaction_handle.theta0;
                    LC::scalar phi0 = M_PI / 180.f * _widget.knot_interaction_handle.phi0;

                    disp[0] = r * sin(theta0) * cos(phi0);
                    disp[1] = r * sin(theta0) * sin(phi0);
                    disp[2] = r * cos(theta0);
                    
                    // Set the background field to uniform
                    clear_heliknotons();

                    // Embed dimer centered at each given position
                    for (auto& pos : _widget.knot_interaction_handle.pos_gui.position_array) {
                        Eigen::Vector3d p0(pos[0], pos[1], pos[2]);
                        Eigen::Vector3d ph1 = p0 + disp;
                        Eigen::Vector3d ph2 = p0 - disp;
                        translations.emplace_back(ph1);
                        translations.emplace_back(ph2);
                        embed_heliknoton(ph1);
                        embed_heliknoton(ph2);
                    }

                }
            }

            for (int i = 0; i < VOX[0]; i++) {
                for (int j = 0; j < VOX[1]; j++) {
                    for (int k = 0; k < VOX[2]; k++) {
                        full_region.push_back(i + VOX[0] * j + VOX[0] * VOX[1] * k);
                    }
                }
            }

            // Set the relax rate
            solver->SetRate(_widget.knot_interaction_handle.relaxRate);

            // 1.5 Relax for specified initial iterations
            if (_widget.knot_interaction_handle.nPrerelax) {
                // Relax
                if (_widget.knot_interaction_handle.useInitialConditions)
                    solver->DomainRelax(_widget.knot_interaction_handle.nPrerelax, full_region, true, true);
                else
                    solver->Relax(_widget.knot_interaction_handle.nPrerelax, true);
            }

            // Create file to log energies for each knot
            std::string fname_energy = _widget.knot_interaction_handle.FileName() + "_energy.bin";
            std::vector<LC::scalar> energy_array;
            std::string fname_positions = _widget.knot_interaction_handle.FileName() + "_positions.bin";
            std::vector<std::vector<Eigen::Vector3d>> positions_array;

            // Need to extract previous iterations
            if (_widget.knot_interaction_handle.nFrameOffset) {
                // Load previous energy data
                {
                    std::streampos fsize = 0;
                    {
                        std::ifstream tempfile(fname_energy, std::ios::binary);

                        fsize = tempfile.tellg();
                        tempfile.seekg(0, std::ios::end);
                        fsize = tempfile.tellg() - fsize;
                        tempfile.close();
                    }

                    int temp_points = fsize / sizeof(LC::scalar);
                    std::vector<LC::scalar> tempdata(temp_points);

                    // Read in data
                    std::ifstream ifile(fname_energy, std::ios::in | std::ios::binary);
                    if (ifile.is_open() && fsize > 0) {
                        ifile.read((char*)&tempdata[0], fsize);
                        ifile.close();
                    }
                    energy_array = tempdata;
                }
                // Load previous position data
                {
                    std::streampos fsize = 0;
                    {
                        std::ifstream tempfile(fname_positions, std::ios::binary);

                        fsize = tempfile.tellg();
                        tempfile.seekg(0, std::ios::end);
                        fsize = tempfile.tellg() - fsize;
                        tempfile.close();
                    }

                    // Read in data
                    std::ifstream ifile(fname_positions, std::ios::in | std::ios::binary);
                    if (ifile.is_open() && fsize > 0) {
                        // Get the number of frames that will be read
                        int pos_frames_detected;
                        ifile.read((char*)&pos_frames_detected, sizeof(int));
                        // For each frame
                        for (int f = 0; f < pos_frames_detected; f++) {
                            // Read the number of positions
                            int nPos;
                            ifile.read((char*)&nPos, sizeof(int));
                            std::vector<Eigen::Vector3d> tmppos(nPos);
                            ifile.read((char*)&tmppos[0], nPos * sizeof(Eigen::Vector3d));
                            positions_array.push_back(tmppos);
                        }
                        ifile.close();
                    }
                }

            }
            // Save info file about the simulation
            else {
                std::string info_txt = _widget.knot_interaction_handle.InfoFile();
                std::ofstream info_out(info_txt);
                if (info_out.is_open()) {

                    info_out << "CELL = " << std::to_string(CELL[0]) + "X" + std::to_string(CELL[1]) + "X" + std::to_string(CELL[2]) << std::endl;
                    info_out << "Node density = " + std::to_string(VOX[0]/CELL[0]) +
                        "X" + std::to_string(VOX[1] / CELL[1]) +
                        "X" + std::to_string(VOX[2] / CELL[2])
                        << std::endl;
                    info_out << "Vortex line Upsampling = " + std::to_string(_widget.knot_interaction_handle.upsampling_multiplier) + "x" << std::endl;
                    info_out << "Iterations/Frame = " + std::to_string(_widget.knot_interaction_handle.nRelax_per_frame) << std::endl;
                    info_out << "Iterations (Prerelax) = " + std::to_string(_widget.knot_interaction_handle.nPrerelax) << std::endl;

                    for (int i = 0; i < _widget.knot_interaction_handle.switching_data.size(); i++) {
                        float efieldval = _widget.knot_interaction_handle.switching_data[i].efield;
                        int dframes = _widget.knot_interaction_handle.switching_data[i].dframes;
                        info_out << "E (" + std::to_string(i) + ") = " + std::to_string(efieldval) << std::endl;
                        info_out << "dFrames (" + std::to_string(i) + ") = " + std::to_string(dframes) << std::endl;
                    }
                    if (_widget.knot_interaction_handle.cycleEfield)
                        info_out << "Efield cycling = ON" << std::endl;
                    else
                        info_out << "Efield cycling = OFF" << std::endl;

                    info_out << "Rate = " + std::to_string(_widget.knot_interaction_handle.relaxRate) << std::endl;
                    info_out << "Frames = " + std::to_string(_widget.knot_interaction_handle.nFrames);
                    info_out.close();
                }
            }

            // Initialize data once
            uint slice_vortex = data->voxels[0] * data->voxels[1];
            uint vol_vortex = data->voxels[2] * slice_vortex;
            std::array<int, 3> vox_interp_vortex; // Can play around with this parameter
            for (int d = 0; d < 3; d++)
                vox_interp_vortex[d] = VOX[d] * _widget.knot_interaction_handle.upsampling_multiplier;
            uint vol_interp_vortex = vox_interp_vortex[0] * vox_interp_vortex[1] * vox_interp_vortex[2];

            std::unique_ptr<float[]> chi_field(new float[3 * vol_vortex]), field_nn(new float[3 * vol_vortex]),
                valid_fieldf(new float[vol_vortex]),
                field_S(new float[vol_vortex]),
                sample_grid(new float[vol_interp_vortex]),
                field_Q(new float[5 * vol_vortex]),
                field_Ssb(new float[vol_vortex]);
            _vortex_line.valid_field = std::unique_ptr<short[]>(new short[vol_vortex]);

            // Initialize outside loop in vortex interactions code
            LC::ExtendedMC::MarchingCubes mc, mc2;
            mc.set_resolution(data->voxels[0], data->voxels[1], data->voxels[2]);
            mc.init_all();
            mc2.set_resolution(vox_interp_vortex[0], vox_interp_vortex[1], vox_interp_vortex[2]);
            mc2.init_all();

            // Update visual data sizes to use POM() function
            initVisuals();

            // Set location to save the pom image
            _widget.savePOM = true;
            _widget.savePOM_loc = _widget.knot_interaction_handle.PathToFile();
            _image_series.SetCount(_widget.knot_interaction_handle.nFrameOffset);

            // 3 Run the routine
            int current_switch_index = 0;
            int switch_frame_count = 0;
            int num_switches = _widget.knot_interaction_handle.switching_data.size();

            for (int f = 0; f < _widget.knot_interaction_handle.nFrames; f++) {

                int f_eff = _widget.knot_interaction_handle.nFrameOffset + f;

                if (num_switches) {
                    auto current_switch = _widget.knot_interaction_handle.switching_data[current_switch_index];

                    // Set the first switch
                    if (f == 0) {
                        solver->SetVoltage(current_switch.efield * data->cell_dims[2], 500);
                    }

                    if (current_switch.dframes == switch_frame_count++) {
                        // Move to the next switch
                        current_switch_index++;

                        // If the switch index exceeds the maximum number...
                        if (_widget.knot_interaction_handle.cycleEfield) {
                            current_switch_index = current_switch_index % num_switches;
                            current_switch = _widget.knot_interaction_handle.switching_data[current_switch_index];
                            solver->SetVoltage(current_switch.efield * data->cell_dims[2], 500);
                        }
                        else if (!_widget.knot_interaction_handle.cycleEfield && current_switch_index < num_switches) {
                            current_switch = _widget.knot_interaction_handle.switching_data[current_switch_index];
                            solver->SetVoltage(current_switch.efield * data->cell_dims[2], 500);
                        }
                        // Reset count
                        switch_frame_count = 0;
                    }
                }
                

                // 3a Relax
                solver->Relax(_widget.knot_interaction_handle.nRelax_per_frame, true);
                std::string fname_full = _widget.knot_interaction_handle.FileName() + std::to_string(f_eff);
                // 3b Generate isosurface and save isosurface as ply file
                auto preimage_components = findVortexKnot(true, fname_full, chi_field, field_nn, valid_fieldf, field_S, sample_grid,field_Q, field_Ssb, mc, mc2, vox_interp_vortex);
                energy_array.push_back(solver->TotalEnergy());
                // 3c Find the vortex knots
                // Set mc field
                mc.set_ext_data(valid_fieldf.get());
                auto rvecs = detect_heliknoton_COM(preimage_components,
                    translations);
                if (rvecs.size()) positions_array.push_back(rvecs);
                else LC_WARN("No heliknotons were detected");
                // 3d Generate the POM image
                POM("POM");
                // Update the energy file each iteration
                {
                    std::ofstream ofile(fname_energy, std::ios::out | std::ios::binary);
                    if (ofile.is_open()) {
                        std::size_t nBytes = sizeof(LC::scalar) * energy_array.size();
                        ofile.write((char*)&energy_array[0], nBytes);
                        ofile.close();
                    }
                }

                // Update the position file each iteration
                {
                    std::ofstream ofile(fname_positions, std::ios::out | std::ios::binary);

                    if (ofile.is_open()) {
                        int nf = positions_array.size();
                        if (nf) {
                            ofile.write((char*)&nf, sizeof(int));
                            // For each frame...
                            for (const auto & h_detected : positions_array) {
                                int com_sz = h_detected.size();
                                ofile.write((char*)&com_sz, sizeof(int));
                                ofile.write((char*)&h_detected[0], com_sz * sizeof(Eigen::Vector3d));
                            }
                        }
                        
                        ofile.close();
                    }
                }

                // Check if backup frame and save
                if (f == 0 || (f+1) % _widget.knot_interaction_handle.nBackup_rate == 0) {
                    std::string res = _widget.knot_interaction_handle.FileName() + std::to_string(f_eff) + ".lmt";
                    if (!res.empty()) {
                        _header.writeFile = res;
                        _solver->Export(_header);
                        LC_CORE_INFO("Saved to file {0}", res.c_str());
                        _header.readFile = res;
                    }
                }
            }

            _widget.updateImage = true;
            _widget.savePOM = false;
        }
    }
}

void Sandbox::handleVortexKnotWindow() {
    if (!_widget.showVortexKnotSettings || !_solver) return;

    auto* data = static_cast<Dataset*>(_solver->GetDataPtr());
    if (!data) return;

    ImGui::SetNextWindowSize({ 410.0f, 620.0f }, ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Knot Search Tool", &_widget.showVortexKnotSettings)) {
        ImGui::End();
        return;
    }

    if (ImGui::CollapsingHeader("Search & Tracing Parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushItemWidth(130.0f);

        if (ImGui::InputInt("Max Vortex Components", &_widget.maxVortexComponents)) {
            _widget.maxVortexComponents = (std::max)(1, _widget.maxVortexComponents);
        }
        if (ImGui::InputInt("Min Component Size", &_widget.minComponentSize)) {
            _widget.minComponentSize = (std::max)(1, _widget.minComponentSize);
        }
        if (ImGui::InputInt("Max Vortex Points", &_vortex_line.maxPoints)) {
            _vortex_line.maxPoints = (std::max)(3, _vortex_line.maxPoints);
        }
        if (ImGui::InputInt("Max Conic Angle (deg)", &_vortex_line.maxConic)) {
            _vortex_line.maxConic = std::clamp(_vortex_line.maxConic, 1, 180);
        }

        int knn = _vortex_line.knn;
        if (ImGui::InputInt("Nearest Neighbors (k > 3)", &knn)) {
            _vortex_line.knn = (std::max)(4, knn);
        }

        if (ImGui::InputFloat("Point Separation", &_vortex_line.pointSeparation, 0.05f, 0.25f, "%.3f")) {
            _vortex_line.pointSeparation = (std::max)(1e-4f, _vortex_line.pointSeparation);
        }
        if (ImGui::InputFloat("Completion Threshold", &_widget.knotCompletionDist, 0.05f, 0.25f, "%.3f")) {
            _widget.knotCompletionDist = (std::max)(0.0f, _widget.knotCompletionDist);
        }

        ImGui::PopItemWidth();
    }

    if (ImGui::CollapsingHeader("Refinement & Mesh", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushItemWidth(130.0f);

        if (ImGui::InputInt("Initial Cutoff Points", &_widget.knotInitialCutoffIterations)) {
            _widget.knotInitialCutoffIterations = (std::max)(0, _widget.knotInitialCutoffIterations);
        }
        if (ImGui::InputInt("Refinement Iterations", &_widget.knotRefinementIterations)) {
            _widget.knotRefinementIterations = (std::max)(0, _widget.knotRefinementIterations);
        }
        if (ImGui::InputInt("Post-Process Upsampling", &_vortex_line.upSample)) {
            _vortex_line.upSample = (std::max)(1, _vortex_line.upSample);
        }
        if (ImGui::InputFloat("Tube Radius", &_vortex_line.tubeRadius, 0.05f, 0.25f, "%.3f")) {
            _vortex_line.tubeRadius = (std::max)(1e-4f, _vortex_line.tubeRadius);
        }

        ImGui::PopItemWidth();

        if (ImGui::Button("Find Components", ImVec2(150.0f, 0.0f))) {
            findVortexKnotComponents();
        }
        ImGui::SameLine();
        ImGui::Checkbox("Draw Surfaces", &_widget.drawSurfaces);

        if (_widget.generateKnots) {
            generateIsosurface();
        }
    }

    const ImGuiTableFlags knotTableFlags = ImGuiTableFlags_BordersInnerV |
        ImGuiTableFlags_RowBg |
        ImGuiTableFlags_SizingStretchProp;

    if (!_vortexKnot.empty() && ImGui::CollapsingHeader("Detected Knots", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::BeginTable("RawKnotsTable", 4, knotTableFlags)) {
            ImGui::TableSetupColumn("Knot", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("Visible", ImGuiTableColumnFlags_WidthFixed, 55.0f);
            ImGui::TableSetupColumn("Color", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableSetupColumn("##Remove", ImGuiTableColumnFlags_WidthFixed, 32.0f);
            ImGui::TableHeadersRow();

            int it = 0;
            for (auto iter = _vortexKnot.begin(); iter != _vortexKnot.end(); ++iter, ++it) {
                ImGui::PushID(it);
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                ImGui::Text("Knot K%d%s", it, iter->ribbon ? " (Ribbon)" : "");

                ImGui::TableSetColumnIndex(1);
                ImGui::Checkbox("##draw", &iter->draw);

                ImGui::TableSetColumnIndex(2);
                if (ImGui::ColorEdit4("##color", &iter->knotColor[0], ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
                    iter->UpdateColor();
                }

                ImGui::TableSetColumnIndex(3);
                if (ImGui::SmallButton("X")) {
                    iter->draw = false;
                    _vortexKnot.erase(iter);
                    if (static_cast<std::size_t>(it) < _vortex_line.components.size()) {
                        _vortex_line.components.erase(_vortex_line.components.begin() + it);
                    }
                    ImGui::PopID();
                    break;
                }

                ImGui::PopID();
            }
            ImGui::EndTable();
        }
    }

    if (!_processedVortexKnot.empty() && ImGui::CollapsingHeader("Processed Knots", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::BeginTable("ProcessedKnotsTable", 4, knotTableFlags)) {
            ImGui::TableSetupColumn("Knot", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("Visible", ImGuiTableColumnFlags_WidthFixed, 55.0f);
            ImGui::TableSetupColumn("Color", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableSetupColumn("##Remove", ImGuiTableColumnFlags_WidthFixed, 32.0f);
            ImGui::TableHeadersRow();

            int it = 0;
            for (auto iter = _processedVortexKnot.begin(); iter != _processedVortexKnot.end(); ++iter, ++it) {
                ImGui::PushID(1000 + it);
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                ImGui::Text("Processed K%d", it);

                ImGui::TableSetColumnIndex(1);
                ImGui::Checkbox("##draw_proc", &iter->draw);

                ImGui::TableSetColumnIndex(2);
                if (ImGui::ColorEdit4("##color_proc", &iter->knotColor[0], ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
                    iter->UpdateColor();
                }

                ImGui::TableSetColumnIndex(3);
                if (ImGui::SmallButton("X")) {
                    iter->draw = false;
                    _processedVortexKnot.erase(iter);
                    ImGui::PopID();
                    break;
                }

                ImGui::PopID();
            }
            ImGui::EndTable();
        }
    }

    if (ImGui::CollapsingHeader("Export Knots", ImGuiTreeNodeFlags_DefaultOpen)) {
        constexpr std::size_t bufSize = 128;
        static char exportDirBuf[bufSize] = "../data";
        static char exportNameBuf[bufSize] = "";
        static float knotPointDensity = 1.0f;

        ImGui::InputText("Directory##Knot", exportDirBuf, bufSize);
        ImGui::InputText("Base Name##Knot", exportNameBuf, bufSize);

        ImGui::SetNextItemWidth(110.0f);
        if (ImGui::InputFloat("Sample Density", &knotPointDensity, 0.1f, 0.5f, "%.2f")) {
            knotPointDensity = (std::max)(0.05f, knotPointDensity);
        }

        if (ImGui::Button("Export Knots", ImVec2(140.0f, 0.0f)) && data->directors) {
            const std::filesystem::path exportDir(exportDirBuf);
            const std::string baseName = std::string(exportNameBuf).empty() ? "knot" : std::string(exportNameBuf);

            const auto& cellf = data->cell_dims;
            const int nx = data->voxels[0];
            const int ny = data->voxels[1];
            const int nz = data->voxels[2];
            const std::size_t vol = static_cast<std::size_t>(nx) * ny * nz;

            if (nx > 1 && ny > 1 && nz > 1) {
                const float dr = static_cast<float>(cellf[0]) / static_cast<float>(nx - 1);
                const LC::scalar minPointDist = 2.0 * knotPointDensity * std::sqrt(3.0) * dr;
                const std::string cumulantBase = (exportDir / (baseName + "_K_vortex.bin")).generic_string();

                bool firstExported = true;
                int it = 0;

                for (auto Ki = _vortexKnot.begin(); Ki != _vortexKnot.end(); ++Ki, ++it) {
                    if (Ki->ribbon) continue;

                    const std::string objPrefix = (exportDir / (baseName + "_K" + std::to_string(it))).generic_string();

                    const std::size_t Nv = Ki->surface.vertices.size();
                    const std::size_t Ni = Ki->surface.indices.size() / 3;
                    std::vector<MeshLib::PNCVertex<float>> Kiv(Nv);
                    std::vector<MeshLib::Triangle> Kitrigs(Ni);

                    for (std::size_t i = 0; i < Nv; ++i) {
                        const auto& Kip = Ki->surface.vertices[i].position;
                        const auto& Kin = Ki->surface.vertices[i].normal;
                        const auto& Kic = Ki->surface.vertices[i].color;
                        Kiv[i].position = { Kip[0], Kip[1], Kip[2] };
                        Kiv[i].normal = { Kin[0], Kin[1], Kin[2] };
                        Kiv[i].color = { Kic[0], Kic[1], Kic[2], Kic[3] };
                    }

                    for (std::size_t i = 0; i < Ni; ++i) {
                        Kitrigs[i][0] = Ki->surface.indices[3 * i];
                        Kitrigs[i][1] = Ki->surface.indices[3 * i + 1];
                        Kitrigs[i][2] = Ki->surface.indices[3 * i + 2];
                    }

                    auto exportedData = exportVortexKnot(Kiv, Kitrigs, minPointDist, objPrefix + "_vortex.bin", false);
                    exportVortexKnot(Kiv, Kitrigs, minPointDist, cumulantBase, !firstExported);
                    firstExported = false;

                    const auto& pdata = std::get<0>(exportedData);
                    if (pdata.empty()) continue;

                    std::vector<Eigen::Vector2d> odata(pdata.size());

                    for (std::size_t p = 0; p < pdata.size(); ++p) {
                        const int ix = std::clamp(static_cast<int>(std::round((pdata[p].x() / cellf[0] + 0.5) * (nx - 1))), 0, nx - 1);
                        const int iy = std::clamp(static_cast<int>(std::round((pdata[p].y() / cellf[1] + 0.5) * (ny - 1))), 0, ny - 1);
                        const int iz = std::clamp(static_cast<int>(std::round((pdata[p].z() / cellf[2] + 0.5) * (nz - 1))), 0, nz - 1);

                        const std::size_t id = static_cast<std::size_t>(ix)
                            + static_cast<std::size_t>(iy) * nx
                            + static_cast<std::size_t>(iz) * nx * ny;

                        const LC::scalar nxVal = data->directors[id];
                        const LC::scalar nyVal = data->directors[id + vol];
                        const LC::scalar nzVal = std::clamp(data->directors[id + 2 * vol], LC::scalar(-1.0), LC::scalar(1.0));

                        odata[p][0] = std::acos(nzVal);
                        odata[p][1] = std::atan2(nyVal, nxVal);
                    }

                    std::ofstream ofile(objPrefix + "_directors.bin", std::ios::out | std::ios::binary);
                    if (ofile.is_open()) {
                        const std::int32_t osz = static_cast<std::int32_t>(odata.size());
                        ofile.write(reinterpret_cast<const char*>(&osz), sizeof(std::int32_t));
                        ofile.write(reinterpret_cast<const char*>(odata.data()), osz * sizeof(Eigen::Vector2d));
                    }
                }
            }
        }
    }

    ImGui::End();
}



void Sandbox::nematicPreimageManager() {

    std::vector <std::array<float, 2>> remove_list;

    for (auto& p : _nematicPreimages)
        if (!p.opentab) {
            // Set to not draw (bad things will happen otherwise):
            p.draw = false;
            remove_list.push_back({ p.theta, p.phi });
        }

    // Now we can remove
    for (const auto& p : remove_list) {
        _nematicPreimages.remove(S2Fiber{ p[0], p[1] });

    }

    // Check if using global alpha
    if (_widget.global_preimage_alpha)
        for (auto& p : _nematicPreimages)
            p.alpha = _widget.preimage_alpha;

    if (_nematicPreimages.size())
        ImGui::Text("(Theta, Phi)");
    char txtbuffer[30];

    if (ImGui::BeginTabBar("Selected preimages", ImGuiTabBarFlags_Reorderable | ImGuiTabBarFlags_AutoSelectNewTabs | ImGuiTabBarFlags_FittingPolicyScroll)) {
        for (auto& p : _nematicPreimages) {

            int n2 = sprintf(txtbuffer, "(%.0f,%.0f)", p.theta, p.phi);
            // Copy to string
            std::string txt;
            for (int i = 0; i < n2; i++) txt += txtbuffer[i];

            if (p.opentab && ImGui::BeginTabItem(txt.c_str(), &p.opentab, ImGuiTabItemFlags_None))
            {
                ImGui::Text("Preimage %s", txt.c_str());
                // Allow isoLevel and alpha to be adjusted
                ImGui::SliderFloat(("Isovalue ##" + txt).c_str(), &p.isoLevel, 0.0f, 0.3f);
                if (!_widget.global_preimage_alpha)
                    ImGui::SliderFloat(("Iso Alpha ##" + txt).c_str(), &p.alpha, 0.0f, 1.0f);

                ImGui::InputInt(("Mesh subdivisions##" + txt).c_str(), &p.subdivisions);
                ImGui::Checkbox("Pontryagin construction", &p.pontryagin);
                ImGui::SameLine();
                ImGui::Checkbox("Invert normals", &p.normal_inversion);
                ImGui::SameLine();
                ImGui::Checkbox("Cull faces", &p.cullFaces);
                ImGui::SameLine();
                ImGui::Checkbox("Domain inversion", &p.domain_inversion);
                ImGui::Checkbox("Draw surface", &p.draw);

                ImGui::EndTabItem();
            }
        }
        ImGui::EndTabBar();
    }

}

void Sandbox::handleZProfileWindow() {

    if (_widget.showZProfileWindow) {

        ImGui::Begin("Z-profile", &_widget.showZProfileWindow);

        Dataset* data = (Dataset*)(_solver->GetDataPtr());

        _widget.generateProfile = ImGui::Button("Generate z-profile");
        ImGui::SameLine();
        ImGui::Checkbox("Draw profile", &_widget.drawProfile);
        
        ImGui::SliderInt("Chain length", &_widget.chain_units, 1, 100);
        ImGui::SliderFloat("Omega bar", &_widget.omegabar, 0.0f, 1.0f);

        ImGui::SliderFloat("Thermal randomization", &_widget.thermalRandomization, 0.0f, M_PI / 2.0f);
        ImGui::SliderFloat("Expansion coefficient", &_widget.expansionCoeff, 1.0f, 3.0f);
        ImGui::SliderFloat("Contraction coefficient", &_widget.contractionCoeff, 0.5f, 1.0f);
        ImGui::SliderFloat("Temperature (K)", &_widget.temperature, 250.0f, 450.0f);
        ImGui::SliderFloat("Inversion Temperature (K)", &_widget.inversion_temp, 250.0f, 400.0f);

        if (_widget.generateProfile) {
            _zprofile = Zprofile{};
            _zprofile->expansion_factor = _widget.expansionCoeff;
            _zprofile->contraction_factor = _widget.contractionCoeff;
            _zprofile->GenerateProfile(data->directors.get(), data->voxels, data->cell_dims,
                _widget.thermalRandomization, _widget.chain_units, _widget.temperature,
            _widget.omegabar, _widget.inversion_temp);

            // For each director on the plane, update the
            // boundary conditions based on the tangent at that point
            int V = data->voxels[0];
            int U = data->voxels[1];
            int W = data->voxels[2];
            uint32_t vol = V * U * W;
            for (int x = 1; x < V-1; x++) {
                for (int y = 1; y < U-1; y++) {
                    // Compute the local tangent
                    auto y1 = _zprofile->graph->data[x + V * (y-1)];
                    auto y2 = _zprofile->graph->data[x + V * (y+1)];
                    auto x1 = _zprofile->graph->data[(x-1) + V * y];
                    auto x2 = _zprofile->graph->data[(x+1) + V * y];
                    Vector3 u = x2 - x1;
                    Vector3 v = y2 - y1;
                    Vector3 n = Math::cross(u, v);
                    n = n / n.length();
                    // Copy to director data
                    data->directors[x + V * y + V * U * (W - 1)] = n.x();
                    data->directors[x + V * y + V * U * (W - 1) + vol] = n.y();
                    data->directors[x + V * y + V * U * (W - 1) + vol * 2] = n.z();
                }
            }
        }

        if (ImGui::RadioButton("x-line", &_widget.radioZProfileAxis, 0))
            _widget.radioZProfileIndex = -1;
        ImGui::SameLine();
        if (ImGui::RadioButton("y-line", &_widget.radioZProfileAxis, 1))
            _widget.radioZProfileIndex = -1;
        ImGui::SameLine();
        ImGui::PushItemWidth(150.0f);
        if (_widget.radioZProfileIndex == -1)
            _widget.radioZProfileIndex = (data->voxels[_widget.radioZProfileAxis]-1) / 2;
        ImGui::SliderInt("Index", &_widget.radioZProfileIndex, 0, data->voxels[_widget.radioZProfileAxis] - 1);
        ImGui::PopItemWidth();

        if (_zprofile) {
            // Plot of the z-profile along x slice

            int vox[] = { _zprofile->graph->vox_x, _zprofile->graph->vox_x };
            int complementAxis = (_widget.radioZProfileAxis + 1) % 2;

            int vox_x = vox[_widget.radioZProfileAxis];
            int vox_y = vox[complementAxis];

            std::unique_ptr<float[]> ax(new float[vox_x]);
            std::unique_ptr<float[]> line_scan(new float[vox_y]);


            float dx = data->cell_dims[_widget.radioZProfileAxis] / (data->voxels[_widget.radioZProfileAxis] - 1);
            float dy = data->cell_dims[complementAxis] / (data->voxels[complementAxis] - 1);

            for (int i = 0; i < vox_x; i++) {
                ax[i] = dx * i;
                if (_widget.radioZProfileAxis == 0)
                    line_scan[i] = _zprofile->graph->data[i + _widget.radioZProfileIndex * vox[0]].z();
                else if (_widget.radioZProfileAxis == 1)
                    line_scan[i] = _zprofile->graph->data[_widget.radioZProfileIndex + i * vox[0]].z();
            }

            std::string letter;

            if (_widget.radioZProfileAxis == 0) // x-axis
                letter = 'x';
            else
                letter = 'y';

            std::string plt_title = "Z-profile plot (" + letter + "-axis)";

            auto fautofit = ImPlotAxisFlags_::ImPlotAxisFlags_AutoFit;

            if (ImPlot::BeginPlot(plt_title.c_str(),"L / p", "h / p", ImVec2(-1,0), 0, fautofit, fautofit)) {
                ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
                ImPlot::PlotLine(("Z(" + letter + ")").c_str(), ax.get(), line_scan.get(), vox_x, 0, sizeof(float));
                ImPlot::EndPlot();
            }

            // Draw the plane corresponding to the cross section
            
            _zprofile->box.dims[complementAxis] = 0;
            _zprofile->box.dims[_widget.radioZProfileAxis] = data->cell_dims[_widget.radioZProfileAxis];
            _zprofile->box.dims[2] = data->cell_dims[2];

            _zprofile->box.translation[complementAxis] = dy * (_widget.radioZProfileIndex - data->voxels[complementAxis] + 1) + 0.5f * data->cell_dims[complementAxis];
            _zprofile->box.translation[_widget.radioZProfileAxis] = 0.0f;
            _zprofile->box.translation[2] = 0.0f;
            
        }

        ImGui::End();

    }
}

/*
    Duplicate the volume along the specified direction
    i = 0 : x
    i = 1 : y
    i = 2 : z
*/
void Sandbox::repeatVolume(int i_) {

    auto data = (Dataset*)_solver->GetDataPtr();
    auto solver = (FOFDSolver*)_solver.get();

    std::string sym;
    if (i_ == 0) sym = 'x';
    else if (i_ == 1) sym = 'y';
    else if (i_ == 2) sym = 'z';
    std::string label = "Repeat (" + sym + ")";

    if (ImGui::Button(label.c_str())) {
        auto voxi = data->voxels;
        auto vNew = voxi;
        vNew[i_] *= 2;

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

                        for (int d = 0; d < 3; d++) {
                            if (i_ == 0)
                                nn_new(i + voxi[0] * rep, j, k, d) = nn(i, j, k, d);
                            else if (i_ == 1)
                                nn_new(i, j + voxi[1] * rep, k, d) = nn(i, j, k, d);
                            else if (i_ == 2)
                                nn_new(i, j, k + voxi[2] * rep, d) = nn(i, j, k, d);
                        }

                    }

        LC::scalar v0 = vv(0, 0, voxi[2] - 1);
        data->cell_dims[i_] *= 2;

        // Replace data
        for (int d = 0; d < 3; d++) data->voxels[d] = vNew[d];
        data->directors = std::unique_ptr<LC::scalar[]>(field_nn.release());
        data->voltage = std::unique_ptr<LC::scalar[]>(field_vv.release());
        data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);
        for (int d = 0; d < 3; d++)
            _widget.celldims[d] = data->cell_dims[d];
        solver->SetVoltage(v0);
        initVisuals();
        _widget.updateImageFromLoad = true;
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

        // Only show scan if location has been chosen
        if (_widget.saveNonlin_loc.size()) {
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
        }
        ImGui::End();
    }
}

void Sandbox::handleLCINFOWindow() {
    if (!_widget.showLCINFO || !_solver) return;

    using namespace LC::FrankOseen;

    auto* data = static_cast<Dataset*>(_solver->GetDataPtr());
    if (!data) return;

    ImGui::SetNextWindowSize(ImVec2(480.0f, 520.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("LC Info", &_widget.showLCINFO)) {
        ImGui::End();
        return;
    }

    if (ImGui::CollapsingHeader("Relaxation Convergence", ImGuiTreeNodeFlags_DefaultOpen)) {
        int points = static_cast<int>(_widget.series_x_axis.size());
        if (_widget.updateImageFromLoad) {
            points = _widget.max_graph_points;
        }

        const int safeMaxPoints = static_cast<int>(
            (std::min)(_widget.series_x_axis_vec.size(), _widget.energy_series_vec.size())
            );
        points = std::clamp(points, 0, safeMaxPoints);

        const ImPlotAxisFlags axisFlags = ImPlotAxisFlags_AutoFit;
        if (ImPlot::BeginPlot("Free Energy vs. Iterations", "Iterations", "Free Energy (Kp)", ImVec2(-1.0f, 220.0f), ImPlotFlags_None, axisFlags, axisFlags)) {
            ImPlot::SetNextMarkerStyle(ImPlotMarker_None);
            if (points > 0) {
                ImPlot::PlotLine(
                    "F",
                    _widget.series_x_axis_vec.data(),
                    _widget.energy_series_vec.data(),
                    points,
                    0,
                    sizeof(LC::precision_scalar)
                );
            }
            ImPlot::EndPlot();
        }

        static int maxGraphPointsInput = _widget.max_graph_points;
        ImGui::SetNextItemWidth(110.0f);
        if (ImGui::InputInt("History Buffer Size##LCINFO", &maxGraphPointsInput)) {
            maxGraphPointsInput = (std::max)(10, maxGraphPointsInput);
        }
        ImGui::SameLine();
        if (ImGui::Button("Apply##LCINFO-PLOT")) {
            _widget.max_graph_points = maxGraphPointsInput;
            _widget.series_x_axis_vec.resize(_widget.max_graph_points, 0.0);
            _widget.energy_series_vec.resize(_widget.max_graph_points, 0.0);
        }
    }

    if (ImGui::CollapsingHeader("System & Material Properties", ImGuiTreeNodeFlags_DefaultOpen)) {
        static const std::map<LC_TYPE, std::string> lcMap = LiquidCrystal::Map();
        static const auto relaxMethodMap = Dataset::RelaxMap();

        auto lcIt = lcMap.find(data->lc_type);
        const char* lcName = (lcIt != lcMap.end()) ? lcIt->second.c_str() : "Custom / Unknown";

        auto relaxIt = relaxMethodMap.find(data->relaxKind);
        const char* relaxName = (relaxIt != relaxMethodMap.end()) ? relaxIt->second.c_str() : "Unknown";

        float dV = 0.0f;
        if (data->voltage && data->voxels[0] > 0 && data->voxels[1] > 0 && data->voxels[2] > 0) {
            FOFDSolver::Tensor3 voltage(data->voltage.get(), data->voxels[0], data->voxels[1], data->voxels[2]);
            dV = static_cast<float>(voltage(0, 0, data->voxels[2] - 1) - voltage(0, 0, 0));
        }

        const ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersInnerH |
            ImGuiTableFlags_BordersOuter |
            ImGuiTableFlags_RowBg |
            ImGuiTableFlags_SizingStretchProp;

        if (ImGui::BeginTable("LCInfoTable", 2, tableFlags)) {
            ImGui::TableSetupColumn("Property", ImGuiTableColumnFlags_WidthFixed, 155.0f);
            ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch);

            auto addRow = [](const char* label, const char* fmt, ...) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::TextDisabled("%s", label);
                ImGui::TableSetColumnIndex(1);
                va_list args;
                va_start(args, fmt);
                ImGui::TextV(fmt, args);
                va_end(args);
                };

            addRow("Liquid Crystal", "%s", lcName);
            addRow("Relaxation Method", "%s", relaxName);
            addRow("Iterations", "%d", data->numIterations);
            addRow("Voltage (dV)", "%.3f V", dV);
            addRow("Grid Voxels", "%d x %d x %d", data->voxels[0], data->voxels[1], data->voxels[2]);
            addRow("Cell Dimensions", "%.2f x %.2f x %.2f", data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);
            addRow("Periodic BCs", "[%d, %d, %d]", data->bc[0], data->bc[1], data->bc[2]);

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::AlignTextToFramePadding();
            ImGui::TextDisabled("Topological Charge");
            ImGui::TableSetColumnIndex(1);
            ImGui::AlignTextToFramePadding();
            ImGui::Text("%.3f", _widget.computed_topological_charge);
            ImGui::SameLine();
            if (ImGui::SmallButton("Compute Q") && data->directors) {
                if (checkRelax()) {
                    _widget.computed_topological_charge = LC::Math::ComputeHopfCharge(data->directors.get(), data->voxels);
                }
            }
            ImGui::EndTable();
        }
    }
    ImGui::End();
}

// Find the new positions after relaxation by computing the COM of the n = +- preimages for each heliknoton
std::vector<Eigen::Vector3d> detect_heliknoton_COM(std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>>&preimages, const std::vector<Eigen::Vector3d>& translations) {
        
    std::vector<Eigen::Vector3d> heliknotons_found, displacement_vectors;
    int pimage_cnt = 0;
    for (auto& preimage_components : preimages) {
        ++pimage_cnt;
        std::vector<uint>& indices = std::get<1>(preimage_components);
        std::vector<LC::Math::IsoVertex>& verts = std::get<0>(preimage_components);
        unsigned int nTriangles = indices.size() / 3;
        uint nVert = verts.size();

        std::map<uint, Face> unvisited_list_local;

        // Note that verts currently refers to verts and not global indices
        // Therefore, convert indices data to global index space through vertex positions
        for (int tri = 0; tri < nTriangles; tri++) {
            uint i1 = indices[3 * tri];
            uint i2 = indices[3 * tri + 1];
            uint i3 = indices[3 * tri + 2];

            // Store in map
            unvisited_list_local.insert({ tri, Face(i1, i2, i3) });
        }

        // 2. Separate by components
        auto components = find_all_components_graph(unvisited_list_local, nVert, 20);
        LC_INFO("Components detected for preimage[{1}] = {0}", components.size(), pimage_cnt);

        // 3. Compute the COM for each component
        if (components.size() > 0) {
            std::vector<Eigen::Vector3d> COMs(components.size());
            for (int i = 0; i < components.size(); i++) {
                COMs[i] = Eigen::Vector3d(0., 0., 0.);
                // Compute the COM
                for (auto ci : components[i]) {
                    COMs[i](0) += verts[ci].position[0] / float(components[i].size());
                    COMs[i](1) += verts[ci].position[1] / float(components[i].size());
                    COMs[i](2) += verts[ci].position[2] / float(components[i].size());
                }
            }
            // Sort the COMs such that [0,1] are the closest points to <translation>
            if (components.size() == 2) {
                // 4. Identify the COM with the closest pair of preimages
                Eigen::Vector3d com = COMs[0] - COMs[1];
                displacement_vectors.emplace_back(com);
            }
            else { // Invalid data point return empty vector
                return heliknotons_found;
            }
        }

    }

    // From the displacement vectors computed, compute the symmetrized or antisymmetrized position (whichever is greater)
    if (displacement_vectors.size() == 2) {
        Eigen::Vector3d sum = 0.5 * (displacement_vectors[0] + displacement_vectors[1]);
        Eigen::Vector3d diff = 0.5 * (displacement_vectors[0] - displacement_vectors[1]);
        if (sum.squaredNorm() >= diff.squaredNorm())
            heliknotons_found.emplace_back(sum);
        else
            heliknotons_found.emplace_back(diff);
    }
    return heliknotons_found;
}

void Sandbox::handleModificationWindow() {
    if (!_widget.showModificationWindow || !_solver) return;

    auto* solver = static_cast<FOFDSolver*>(_solver.get());
    auto* data = static_cast<Dataset*>(_solver->GetDataPtr());
    if (!data) return;

    ImGui::SetNextWindowSize({ 430.0f, 680.0f }, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos({ 0.0f, 450.0f }, ImGuiCond_FirstUseEver);

    if (!ImGui::Begin("LC Modification Window", &_widget.showModificationWindow)) {
        ImGui::End();
        return;
    }

    if (ImGui::CollapsingHeader("File I/O", ImGuiTreeNodeFlags_DefaultOpen)) {
        constexpr std::size_t bufSize = 128;
        static char fileDir[bufSize] = "../data/";
        static char loadFile[bufSize] = "";
        static char saveFile[bufSize] = "";

        ImGui::InputText("Directory", fileDir, bufSize);
        ImGui::InputText("Load File", loadFile, bufSize);
        ImGui::InputText("Save File", saveFile, bufSize);

        if (ImGui::Button("Load##ModificationWindow", ImVec2(120.0f, 0.0f))) {
            std::filesystem::path fullPath = std::filesystem::path(fileDir) / loadFile;
            std::string loadName = fullPath.generic_string();

            std::ifstream checkLoadFile(loadName);
            if (!std::string(loadFile).empty() && checkLoadFile.good()) {
                checkLoadFile.close();
                if (checkRelax()) {
                    _header.readFile = loadName;
                    _solver->Import(_header);
                    initVisuals();
                    _widget.updateImageFromLoad = true;
                    LC_INFO("Loaded file {0}", loadName.c_str());
                }
            }
        }

        ImGui::SameLine();
        if (ImGui::Button("Save##ModificationWindow", ImVec2(120.0f, 0.0f))) {
            if (!std::string(saveFile).empty() && checkRelax()) {
                std::filesystem::path fullPath = std::filesystem::path(fileDir) / saveFile;
                std::string saveName = fullPath.generic_string();

                _header.writeFile = saveName;
                _solver->Export(_header);
                _header.readFile = saveName;
                LC_INFO("Saved to file {0}", saveName.c_str());
            }
        }
    }

    if (ImGui::CollapsingHeader("Cell & Topology Generation", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputFloat3("Cell dims", &_widget.celldims[0], "%.2f");
        ImGui::InputInt3("PBCs (x, y, z)", &_widget.boundaries[0]);

        if (ImGui::BeginTable("CellParamsTable", 3, ImGuiTableFlags_SizingStretchSame)) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::SetNextItemWidth(-40.0f);
            ImGui::InputInt("npp", &_widget.npp);

            ImGui::TableSetColumnIndex(1);
            ImGui::SetNextItemWidth(-30.0f);
            ImGui::InputInt("Q", &_widget.topological_charge);

            ImGui::TableSetColumnIndex(2);
            ImGui::SetNextItemWidth(-65.0f);
            ImGui::InputInt("Chirality", &_widget.chirality);
            ImGui::EndTable();
        }

        static constexpr const char* configNames[] = {
            "Heliknoton",
            "Hopfion",
            "Toron",
            "Twistion (T1B2)",
            "Twistion (T2B2)",
            "CF3",
            "CF1",
            "CF1 (Loop)",
            "Lehmann Cluster Line",
            "PQ Torus Knot (Experimental)"
        };

        ImGui::Combo("Ansatz Type", &_widget.hopfion_type, configNames, IM_ARRAYSIZE(configNames));

        static float phaseOffset = 0.0f;
        static std::array<float, 3> heliLim = { 1.0f, 1.0f, 1.0f };
        static int pqP = 2;
        static int pqQ = 3;

        if (_widget.hopfion_type == 0) {
            ImGui::InputFloat3("Heliknoton Scaling", &heliLim[0], "%.2f");
            ImGui::InputFloat("Phase Offset", &phaseOffset, 0.05f, 0.25f, "%.2f");
        }
        else if (_widget.hopfion_type == 9) {
            ImGui::PushItemWidth(110.0f);
            ImGui::InputInt("P##ModifyPQKnot", &pqP);
            ImGui::SameLine();
            ImGui::InputInt("Q##ModifyPQKnot", &pqQ);
            ImGui::PopItemWidth();
        }

        if (ImGui::Button("Generate New", ImVec2(140.0f, 0.0f))) {
            const int Q = _widget.topological_charge;
            const int npp = (std::max)(1, _widget.npp);
            const std::array<int, 3> V = {
                static_cast<int>(_widget.celldims[0] * npp),
                static_cast<int>(_widget.celldims[1] * npp),
                static_cast<int>(_widget.celldims[2] * npp)
            };

            data->Cell(_widget.celldims[0], _widget.celldims[1], _widget.celldims[2]);

            std::array<LC::scalar, 3> sHeliLim = { heliLim[0], heliLim[1], heliLim[2] };
            Dataset::Config cfg;

            switch (_widget.hopfion_type) {
            case 0: cfg = Dataset::Heliknoton(Q, V, data->cell_dims, Q, sHeliLim, { 0.f, 0.f, 0.f }, phaseOffset); break;
            case 1: cfg = Dataset::Hopfion(Q); break;
            case 2: cfg = Dataset::Toron(data->cell_dims, Q); break;
            case 3: cfg = Dataset::Twistion_T1B2(data->cell_dims); break;
            case 4: cfg = Dataset::Twistion_T2B2(data->cell_dims); break;
            case 5: cfg = Dataset::CF3(data->cell_dims[0]); break;
            case 6: cfg = Dataset::CF1(data->cell_dims); break;
            case 7: cfg = Dataset::CF1_Loop(data->cell_dims); break;
            case 8: cfg = Dataset::LehmannClusterLine(data->cell_dims, V); break;
            case 9: cfg = Dataset::PQTorusKnot(pqP, pqQ, data->cell_dims, V, 0.5); break;
            default: cfg = Dataset::Toron(data->cell_dims, Q); break;
            }

            data->chirality = _widget.chirality;
            data->Voxels(V[0], V[1], V[2])
                .Boundaries(_widget.boundaries[0], _widget.boundaries[1], _widget.boundaries[2]);
            data->Configuration(cfg);
            data->numIterations = 0;

            _widget.energy_series.clear();
            _widget.series_x_axis.clear();
            std::fill(_widget.energy_series_vec.begin(), _widget.energy_series_vec.end(), 0.0);
            std::fill(_widget.series_x_axis_vec.begin(), _widget.series_x_axis_vec.end(), 0.0);

            _solver->Init();
            initVisuals();
            _widget.updateImageFromLoad = true;
        }

        ImGui::SameLine();
        if (ImGui::Button("Apply Cell Metadata", ImVec2(150.0f, 0.0f))) {
            data->chirality = _widget.chirality;
            data->Boundaries(_widget.boundaries[0], _widget.boundaries[1], _widget.boundaries[2])
                .Cell(_widget.celldims[0], _widget.celldims[1], _widget.celldims[2]);
            initVisuals();
            _widget.updateImageFromLoad = true;
        }
    }

    const int maxVox = (std::max)({ data->voxels[0], data->voxels[1], data->voxels[2], 1 });

    if (ImGui::CollapsingHeader("Subvolume Region & Cropping", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::SliderInt3("Start (1-based)", &_widget.shrink_interval_begin[0], 1, maxVox);
        ImGui::SliderInt3("End (1-based)", &_widget.shrink_interval_end[0], 1, maxVox);

        for (int d = 0; d < 3; ++d) {
            const int vMax = (std::max)(1, data->voxels[d]);
            _widget.shrink_interval_begin[d] = std::clamp(_widget.shrink_interval_begin[d], 1, (std::max)(1, vMax - 1));
            _widget.shrink_interval_end[d] = std::clamp(_widget.shrink_interval_end[d], _widget.shrink_interval_begin[d] + 1, vMax);
        }

        ImGui::TextDisabled("Active Region: nn(%d:%d, %d:%d, %d:%d)",
            _widget.shrink_interval_begin[0], _widget.shrink_interval_end[0],
            _widget.shrink_interval_begin[1], _widget.shrink_interval_end[1],
            _widget.shrink_interval_begin[2], _widget.shrink_interval_end[2]
        );

        if (ImGui::Button("Crop to Region")) {
            std::array<int, 3> vNew;
            for (int d = 0; d < 3; ++d) {
                vNew[d] = _widget.shrink_interval_end[d] - _widget.shrink_interval_begin[d] + 1;
                data->cell_dims[d] *= static_cast<LC::scalar>(vNew[d]) / static_cast<LC::scalar>(data->voxels[d]);
            }

            const std::size_t newTotal = static_cast<std::size_t>(vNew[0]) * vNew[1] * vNew[2];
            auto field_nn = std::make_unique<LC::scalar[]>(3 * newTotal);
            auto field_vv = std::make_unique<LC::scalar[]>(newTotal);

            if (!data->voltage) {
                solver->SetVoltage(1.0);
            }

            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
            FOFDSolver::Tensor3 vv(data->voltage.get(), data->voxels[0], data->voxels[1], data->voxels[2]);
            FOFDSolver::Tensor4 nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);
            FOFDSolver::Tensor3 vv_new(field_vv.get(), vNew[0], vNew[1], vNew[2]);

            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; ++i) {
                for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; ++j) {
                    for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; ++k) {
                        const int ip = i - _widget.shrink_interval_begin[0] + 1;
                        const int jp = j - _widget.shrink_interval_begin[1] + 1;
                        const int kp = k - _widget.shrink_interval_begin[2] + 1;

                        for (int d = 0; d < 3; ++d) {
                            nn_new(ip, jp, kp, d) = nn(i, j, k, d);
                        }
                        vv_new(ip, jp, kp) = vv(i, j, k);
                    }
                }
            }

            data->voxels = vNew;
            data->directors = std::move(field_nn);
            data->voltage = std::move(field_vv);
            data->en_density = std::make_unique<LC::scalar[]>(newTotal);

            initVisuals();
            _widget.updateImageFromLoad = true;
        }

        ImGui::SameLine();
        repeatVolume(0);
        ImGui::SameLine();
        repeatVolume(1);
        ImGui::SameLine();
        repeatVolume(2);
    }

    if (ImGui::CollapsingHeader("Field Embedding & Transformations", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::SliderInt3("Insert Pos", &_widget.heliknotonInsertLocation[0], 1, maxVox);
        for (int d = 0; d < 3; ++d) {
            _widget.heliknotonInsertLocation[d] = std::clamp(_widget.heliknotonInsertLocation[d], 1, data->voxels[d]);
        }

        if (ImGui::Button("Insert Heliknoton##ModifyLC")) {
            _widget.updateImageFromLoad = true;
            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

            Eigen::Vector3d translation;
            for (int d = 0; d < 3; ++d) {
                translation[d] = (static_cast<LC::scalar>(_widget.heliknotonInsertLocation[d]) / (data->voxels[d] - 1) - 0.5) * data->cell_dims[d];
            }

            const LC::scalar omega0 = 2.0 * M_PI * data->chirality * translation[2];
            auto heliField = LC::Math::Heliknoton(1, data->cell_dims, 1.0, 1.0, { 0.0, 0.0, 0.0 }, omega0, false);

            for (int i = 0; i < data->voxels[0]; ++i) {
                for (int j = 0; j < data->voxels[1]; ++j) {
                    for (int k = 0; k < data->voxels[2]; ++k) {
                        Eigen::Vector3d pos(
                            data->cell_dims[0] * (static_cast<LC::scalar>(i) / (data->voxels[0] - 1) - 0.5),
                            data->cell_dims[1] * (static_cast<LC::scalar>(j) / (data->voxels[1] - 1) - 0.5),
                            data->cell_dims[2] * (static_cast<LC::scalar>(k) / (data->voxels[2] - 1) - 0.5)
                        );
                        pos -= translation;

                        auto res = heliField(pos[0], pos[1], pos[2]);
                        if (res[0] != 0.0 || res[1] != 0.0 || res[2] != 0.0) {
                            nn(i, j, k, 0) = res[0];
                            nn(i, j, k, 1) = res[1];
                            nn(i, j, k, 2) = res[2];
                        }
                    }
                }
            }
        }
        ImGui::SameLine();
        ImGui::Checkbox("Hide Indicator", &_widget.hideHeliknotonInsertIndicator);

        ImGui::Separator();

        if (ImGui::Button("Helical Field")) {
            _widget.updateImage = true;
            auto helical = LC::Math::Planar(2, 1);
            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

            const LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
            const LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
            const LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);

            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; ++i) {
                for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; ++j) {
                    for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; ++k) {
                        const LC::scalar px = -0.5 * data->cell_dims[0] + i * dx;
                        const LC::scalar py = -0.5 * data->cell_dims[1] + j * dy;
                        const LC::scalar pz = -0.5 * data->cell_dims[2] + k * dz - _widget.helicalLayerOffset;

                        auto n = helical(px, py, pz);
                        for (int d = 0; d < 3; ++d) {
                            nn(i, j, k, d) = n[d];
                        }
                    }
                }
            }
        }
        ImGui::SameLine();
        ImGui::SetNextItemWidth(90.0f);
        ImGui::InputFloat("Layer Offset (p)", &_widget.helicalLayerOffset, 0.05f, 0.25f, "%.2f");

        static bool s_uniformFlip = false;
        if (ImGui::Button("Uniform Field")) {
            _widget.updateImage = true;
            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
            const LC::scalar nzVal = s_uniformFlip ? -1.0 : 1.0;

            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; ++i) {
                for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; ++j) {
                    for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; ++k) {
                        nn(i, j, k, 0) = 0.0;
                        nn(i, j, k, 1) = 0.0;
                        nn(i, j, k, 2) = nzVal;
                    }
                }
            }
        }
        ImGui::SameLine();
        ImGui::Checkbox("Invert (-z)", &s_uniformFlip);
        ImGui::SameLine();

        if (ImGui::Button("Flip Chirality") && data->directors) {
            const int nx = data->voxels[0];
            const int ny = data->voxels[1];
            const int nz = data->voxels[2];
            const std::size_t totalVoxels = static_cast<std::size_t>(nx) * ny * nz;

            FOFDSolver::Tensor4 nn(data->directors.get(), nx, ny, nz, 3);
            for (std::size_t idx = 0; idx < totalVoxels / 2; ++idx) {
                const int i = idx % nx;
                const int j = (idx / nx) % ny;
                const int k = idx / (static_cast<std::size_t>(nx) * ny);

                const int i_f = nx - 1 - i;
                const int j_f = ny - 1 - j;
                const int k_f = nz - 1 - k;

                for (int d = 0; d < 3; ++d) {
                    std::swap(nn(i, j, k, d), nn(i_f, j_f, k_f, d));
                }
            }
            _widget.updateImage = true;
        }

        ImGui::Separator();
        ImGui::Text("Z-Plane Translation:");
        const bool arrowUp = ImGui::ArrowButton("Up", ImGuiDir_Up);
        ImGui::SameLine();
        const bool arrowDown = ImGui::ArrowButton("Down", ImGuiDir_Down);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(90.0f);
        ImGui::InputInt("z-planes", &_widget.translationNumber);
        ImGui::SameLine();
        ImGui::Checkbox("Helical", &_widget.helicalTranslation);

        if (arrowUp || arrowDown) {
            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
            const int nz = data->voxels[2];
            const int startz = arrowUp ? _widget.shrink_interval_end[2] - 1 : _widget.shrink_interval_begin[2] - 1;
            const int endz = arrowUp ? _widget.shrink_interval_begin[2] - 1 : _widget.shrink_interval_end[2] - 1;
            const int translate = arrowDown ? -_widget.translationNumber : _widget.translationNumber;

            auto loop = LC::Utility::create_search_list(LC::Utility::irange_pair(startz, endz));
            const LC::scalar dz = data->cell_dims[2] / (nz - 1);
            const LC::scalar dtheta = 2.0 * M_PI * data->chirality * dz;
            const LC::scalar cosNtheta = std::cos(translate * dtheta);
            const LC::scalar sinNtheta = std::sin(translate * dtheta);

            auto wrapZ = [nz](int kVal) -> int {
                return ((kVal % nz) + nz) % nz;
                };

            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; ++i) {
                for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; ++j) {
                    for (const auto& k : loop) {
                        const int kup = wrapZ(k + translate);
                        const int kdown = wrapZ(k - translate);

                        for (int d = 0; d < 3; ++d) {
                            nn(i, j, kup, d) = nn(i, j, k, d);
                            nn(i, j, k, d) = nn(i, j, kdown, d);
                        }

                        if (_widget.helicalTranslation) {
                            const LC::scalar n1x = nn(i, j, kup, 0);
                            const LC::scalar n1y = nn(i, j, kup, 1);
                            const LC::scalar n2x = nn(i, j, k, 0);
                            const LC::scalar n2y = nn(i, j, k, 1);

                            nn(i, j, kup, 0) = n1x * cosNtheta - n1y * sinNtheta;
                            nn(i, j, kup, 1) = n1x * sinNtheta + n1y * cosNtheta;
                            nn(i, j, k, 0) = n2x * cosNtheta - n2y * sinNtheta;
                            nn(i, j, k, 1) = n2x * sinNtheta + n2y * cosNtheta;
                        }
                    }
                }
            }
            _widget.updateImage = true;
        }

        ImGui::Separator();
        ImGui::SetNextItemWidth(100.0f);
        ImGui::InputInt("Scale Factor", &_widget.interpolate);
        ImGui::SameLine();
        if (ImGui::Button("Interpolate Grid") && _widget.interpolate > 1) {
            const int factor = _widget.interpolate;
            const std::size_t oldTotal = static_cast<std::size_t>(data->voxels[0]) * data->voxels[1] * data->voxels[2];
            const std::size_t newTotal = oldTotal * factor * factor * factor;

            auto field_nn = std::make_unique<LC::scalar[]>(3 * newTotal);
            auto field_vv = std::make_unique<LC::scalar[]>(newTotal);

            if (!data->voltage) {
                solver->SetVoltage(1.0);
            }

            LC::Math::Interp3(data->directors.get(), field_nn.get(), data->voxels, { factor, factor, factor }, 3);
            LC::Math::Interp3(data->voltage.get(), field_vv.get(), data->voxels, { factor, factor, factor }, 1);

            for (int d = 0; d < 3; ++d) {
                data->voxels[d] *= factor;
            }

            data->directors = std::move(field_nn);
            data->voltage = std::move(field_vv);

            for (std::size_t i = 0; i < newTotal; ++i) {
                LC::scalar normSq = 0.0;
                for (int d = 0; d < 3; ++d) {
                    const LC::scalar comp = data->directors[i + d * newTotal];
                    normSq += comp * comp;
                }
                if (normSq > 0.0) {
                    const LC::scalar invNorm = 1.0 / std::sqrt(normSq);
                    for (int d = 0; d < 3; ++d) {
                        data->directors[i + d * newTotal] *= invNorm;
                    }
                }
                else {
                    data->directors[i + 2 * newTotal] = 1.0;
                }
            }

            data->en_density = std::make_unique<LC::scalar[]>(newTotal);
            solver->Normalize();
            initVisuals();
            updateColor();
        }
    }
    ImGui::End();
}

void Sandbox::handleInteractionEnergyWindow() {

    if (!_widget.showInteractionEnergyWindow)
        return;

    ImGui::Begin("Interaction Energy Landscape", &_widget.showInteractionEnergyWindow);

    FOFDSolver* solver = (FOFDSolver*)(_solver.get());
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    // 0 -> spherical interaction energy
    // 1 -> radial interaction energy
    static int radioInteractionType = 0;
    static float a_axis = 1.f;
    static float b_axis = 1.f;
    static float c_axis = 1.f;
    ImGui::Separator();
    ImGui::RadioButton("Spherical", &radioInteractionType, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Volumetric", &radioInteractionType, 2);
    ImGui::SameLine();
    ImGui::RadioButton("Radial", &radioInteractionType, 1);
    ImGui::PushItemWidth(100.f);
    ImGui::InputFloat("x radius", &a_axis);
    ImGui::SameLine();
    ImGui::InputFloat("y radius", &b_axis);
    ImGui::SameLine();
    ImGui::InputFloat("z radius", &c_axis);
    ImGui::PopItemWidth();
    ImGui::Separator();

    // Radial parameters

    static float start_dist = 2.f;
    static float end_dist = 5.f;
    static int r_points = 10;
    static float phi0 = 0.f;
    static float theta0 = 0.;
    static int Navg = 1;
    static float EField = 0.5f;
    static float voltage = EField * _widget.interactionCellSize[2];


    ImGui::PushItemWidth(100.f);
    if (radioInteractionType == 0) {
        ImGui::InputFloat("Separation Distance", &_widget.separationDistance);
        ImGui::Checkbox("Use symmetries", &_widget.interactionSymmetry);
        ImGui::InputInt("Theta points", &_widget.interactionThetaPoints);
        ImGui::SameLine();
        ImGui::InputInt("Phi points", &_widget.interactionPhiPoints);
    }
    else if (radioInteractionType == 1) {
        ImGui::InputFloat("Initial r", &start_dist);
        ImGui::InputFloat("End r", &end_dist);
        ImGui::InputInt("r points", &r_points);
        // Make theta and phi interchangable in the future
        ImGui::InputFloat("phi0 (deg)", &phi0);
        ImGui::InputFloat("theta0 (deg)", &theta0);
        ImGui::InputFloat("Interaction offset (deg)", &_widget.interactionOmegaOffset);
    }
    else if (radioInteractionType == 2) {
        ImGui::InputFloat("Initial r", &start_dist);
        ImGui::InputFloat("End r", &end_dist);
        ImGui::InputInt("r points", &r_points);
        ImGui::Checkbox("Use symmetries", &_widget.interactionSymmetry);
        ImGui::InputInt("Theta points", &_widget.interactionThetaPoints);
        ImGui::SameLine();
        ImGui::InputInt("Phi points", &_widget.interactionPhiPoints);
    }
    ImGui::SameLine();
    ImGui::InputInt("Node density (per pitch)", &_widget.interactionNPP);
    ImGui::Checkbox("Single Heliknoton", &_widget.singleInteractionHeliknoton);
    ImGui::PushItemWidth(150.f);
    ImGui::InputFloat3("Cell size", &_widget.interactionCellSize[0]);
    ImGui::PopItemWidth();
    if (ImGui::InputFloat("E-field", &EField)) {
        voltage = EField * _widget.interactionCellSize[2];
    }
    // Create a voltage version as well
    ImGui::SameLine();
    if (ImGui::InputFloat("Voltage", &voltage)) {
        EField = voltage / _widget.interactionCellSize[2];
    }
    ImGui::InputInt("Relaxation iterations", &_widget.interactionIterations);
    ImGui::PopItemWidth();

    // E-field
    const size_t interaction_path_buffer_size = 128;
    static char interaction_path[interaction_path_buffer_size] = "../data/";

    ImGui::InputText("Interaction folder", interaction_path, interaction_path_buffer_size);

    const size_t interaction_fname_buffer_size = 128;
    static char interaction_fname[interaction_fname_buffer_size] = "interaction.bin";

    ImGui::InputText("Interaction file name", interaction_fname, interaction_fname_buffer_size);


    //std::string saveName = "D:/dev/lclab2/data/interactions/" + std::string(interaction_fname);
    std::string saveName = interaction_path + std::string(interaction_fname);

    if (ImGui::Button("Heliknoton interaction landscape")) {

        // Get all points belonging to the heliknoton
        // ============================================
        unsigned int vol = data->voxels[0] * data->voxels[1] * data->voxels[2];
        unsigned int bulk_vol = (data->voxels[0] - 2) * (data->voxels[1] - 2) * (data->voxels[2] - 2);
        typedef Eigen::Vector3d Position;
        typedef Eigen::Vector3d Director;

        struct Point {
            Point() = default;
            Point(const Position& p, const Director& d) : position(p), director(d) {}
            Position position; Director director;
        };

        std::vector<Point> heliknoton;
        heliknoton.reserve(vol);

        Eigen::Matrix3d handedness0;
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                handedness0(a, b) = 0.0;
            }
        }

        // Helical background function of the cell
        auto helical = LC::Math::Planar(2, 1);

        // Normalized handedness to -1
        handedness0(2, 2) = -1.0;

        Position center_of_mass(0., 0., 0);

        // Extract the heliknoton
        for (int i = 0; i < data->voxels[0]; i++) {
            for (int j = 0; j < data->voxels[1]; j++) {
                for (int k = 0; k < data->voxels[2]; k++) {
                    // Diffmag for helical axis
                    unsigned int idx = i + data->voxels[0] * j + data->voxels[0] * data->voxels[1] * k;

                    Position pos;
                    pos[0] = data->cell_dims[0] * ((LC::scalar)i / (data->voxels[0] - 1) - 0.5);
                    pos[1] = data->cell_dims[1] * ((LC::scalar)j / (data->voxels[1] - 1) - 0.5);
                    pos[2] = data->cell_dims[2] * ((LC::scalar)k / (data->voxels[2] - 1) - 0.5);

                    Director dir;
                    dir[0] = data->directors[idx];
                    dir[1] = data->directors[idx + vol];
                    dir[2] = data->directors[idx + 2 * vol];

                    //auto dir_bg_arr = helical(pos[0], pos[1], pos[2]);
                    //Director dir_bg(dir_bg_arr[0], dir_bg_arr[1], dir_bg_arr[2]);

                    // Compute dot product
                    //LC::scalar absprodsq = abs(dir_bg.dot(dir));

                    // Get whatever is inside the ellipsoid
                    if (pow(pos[0] / a_axis, 2) + pow(pos[1] / b_axis, 2) + pow(pos[2] / c_axis, 2) <= 1)
                        heliknoton.push_back({ pos, dir });

                }
            }
        }

        if (!heliknoton.empty())
            center_of_mass = center_of_mass / heliknoton.size();
        else {
            LC_ERROR("Empty heliknoton!!!");
            return;
        }
        // Translate the heliknoton center of mass to the origin
        for (auto& p : heliknoton) {
            p.position = p.position - center_of_mass;
        }

        // Get a bounding box on the heliknoton size
        Position lower_bound(0., 0., 0.), upper_bound(0., 0., 0.);
        for (auto& p : heliknoton) {
            for (int i = 0; i < 3; i++) {
                if (lower_bound[i] > p.position[i])
                    lower_bound[i] = p.position[i];
                if (upper_bound[i] < p.position[i])
                    upper_bound[i] = p.position[i];
            }
        }

        // Pad the background box size
        const LC::scalar background_padding = 1.15;

        // ============================================

        // Create a new volume
        std::array<LC::scalar, 3> cdims;
        if (_widget.forceInteractionCellSize) {
            for (int i = 0; i < 3; i++)
                cdims[i] = _widget.interactionCellSize[i];
        }
        else {
            cdims[0] = 2. * background_padding * (upper_bound[0] - lower_bound[0]) + 0.5 * _widget.separationDistance;
            cdims[1] = 2. * background_padding * (upper_bound[1] - lower_bound[1]) + 0.5 * _widget.separationDistance;
            cdims[2] = 2. * background_padding * (upper_bound[2] - lower_bound[2]) + 0.5 * _widget.separationDistance;
        }

        LC_INFO("Cell dims -> ({0},{1},{2})", cdims[0], cdims[1], cdims[2]);

        std::array<int, 3> vNew = { int(_widget.interactionNPP * cdims[0]), int(_widget.interactionNPP * cdims[1]), int(_widget.interactionNPP * cdims[2]) };
        data->directors = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * vNew[0] * vNew[1] * vNew[2]]);
        std::unique_ptr<LC::scalar[]> nn_data_temp(new LC::scalar[3 * vNew[0] * vNew[1] * vNew[2]]);
        data->voltage = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);
        data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);

        FOFDSolver::Tensor4 nn(data->directors.get(), vNew[0], vNew[1], vNew[2], 3);
        FOFDSolver::Tensor4 nn_temp(nn_data_temp.get(), vNew[0], vNew[1], vNew[2], 3);
        // Change voxels and cell dims
        data->voxels = vNew;
        data->cell_dims = cdims;

        LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
        LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
        LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);

        auto clear_heliknotons = [&]() {
            // Initialize to the uniform helical bg
            for (int i = 0; i < data->voxels[0]; i++) {
                for (int j = 0; j < data->voxels[1]; j++) {
                    for (int k = 0; k < data->voxels[2]; k++) {
                        // New zero point due to translation of center of mass
                        LC::scalar z = -0.5 * data->cell_dims[2] + k * dz;
                        auto bg = helical(0., 0., z);

                        for (int d = 0; d < 3; d++) {
                            nn(i, j, k, d) = bg[d];
                        }
                    }
                }
            }
        };

        // Choose to relax from a relaxed heliknoton (1) or the heliknoton ansatz (0)
#define RELAXED_HELIKNOTON 0

        auto embed_heliknoton = [&](Eigen::Vector3d translation_vector) {

#if RELAXED_HELIKNOTON

            std::set<uint32_t> update_list;

            // reset nn temp
            for (int i = 0; i < vNew[0]; i++)
                for (int j = 0; j < vNew[1]; j++)
                    for (int k = 0; k < vNew[2]; k++)
                        for (int d = 0; d < 3; d++)
                            nn_temp(i, j, k, d) = 0.;



            // Embed the heliknoton in the new volume
            for (const auto& p : heliknoton) {

                // Apply a rotation + translation to heliknoton position
                LC::scalar theta_z = 2. * M_PI * data->chirality * translation_vector[2];
                LC::scalar ctheta = cos(theta_z);
                LC::scalar stheta = sin(theta_z);

                Position newPos;
                Director newDir;

                // While the heliknoton is centered at the origin,
                // apply the z-translation effect to the in plane director positions and total translation
                newPos[0] = ctheta * p.position[0] - stheta * p.position[1] + translation_vector[0];
                newPos[1] = stheta * p.position[0] + ctheta * p.position[1] + translation_vector[1];
                //newPos[0] = p.position[0] + translation_vector[0];
                //newPos[1] = p.position[1] + translation_vector[1];
                newPos[2] = p.position[2] + translation_vector[2];

                // Rotate director due to z translation as well
                newDir[0] = ctheta * p.director[0] - stheta * p.director[1];
                newDir[1] = stheta * p.director[0] + ctheta * p.director[1];
                newDir[2] = p.director[2];

                // Transform new position to index space
                int i = (newPos[0] + 0.5 * data->cell_dims[0]) / dx;
                int j = (newPos[1] + 0.5 * data->cell_dims[1]) / dy;
                int k = (newPos[2] + 0.5 * data->cell_dims[2]) / dz;

                // Add index to update list
                update_list.insert(i + vNew[0] * j + vNew[0] * vNew[1] * k);

                // add translated + rotated point
                for (int d = 0; d < 3; d++)
                    nn_temp(i, j, k, d) += newDir[d];
            }

            for (const auto& idx : update_list) {

                // Transform new position to index space
                int k = idx / (vNew[0] * vNew[1]);
                int j = (idx - k * vNew[0] * vNew[1]) / vNew[0];
                int i = idx - j * vNew[0] - k * vNew[0] * vNew[1];

                // Normalize the points
                LC::scalar nn_length = 0.0;
                for (int d = 0; d < 3; d++)
                    nn_length += nn_temp(i, j, k, d) * nn_temp(i, j, k, d);
                nn_length = sqrt(nn_length);

                for (int d = 0; d < 3; d++)
                    nn(i, j, k, d) = nn_temp(i, j, k, d) / nn_length;
            }

            // Turn update list into a vector
            std::vector<uint32_t> update_list_vec(update_list.size());
            uint32_t count = 0;
            for (auto idx : update_list)
                update_list_vec[count++] = idx;

            return update_list_vec;
#else // Initialize from initial conditions

            LC::scalar omega0 = 2. * M_PI * data->chirality * translation_vector[2] + _widget.interactionOmegaOffset * M_PI / 180.f;
            auto heli_vfield = LC::Math::Heliknoton(1, data->cell_dims, 1., 1., { 0.,0.,0. }, omega0, false);

            for (int i = 0; i < vNew[0]; i++)
                for (int j = 0; j < vNew[1]; j++)
                    for (int k = 0; k < vNew[2]; k++) {
                        // Current coords
                        Position pos;
                        pos[0] = data->cell_dims[0] * ((LC::scalar)i / (data->voxels[0] - 1) - 0.5);
                        pos[1] = data->cell_dims[1] * ((LC::scalar)j / (data->voxels[1] - 1) - 0.5);
                        pos[2] = data->cell_dims[2] * ((LC::scalar)k / (data->voxels[2] - 1) - 0.5);

                        pos = pos - translation_vector;

                        auto res = heli_vfield(pos[0], pos[1], pos[2]);
                        if (res[0] != 0. || res[1] != 0. || res[2] != 0.) {
                            nn(i, j, k, 0) = res[0];
                            nn(i, j, k, 1) = res[1];
                            nn(i, j, k, 2) = res[2];
                        }
                    }
#endif

        };

        std::vector<Eigen::Vector3d> translations;

        // Spherical interaction landscape
        struct InteractionPotential {
            LC::scalar theta;
            LC::scalar phi;
            LC::scalar separation;
            LC::scalar energy;
            LC::scalar en_HH;
            LC::scalar en_H;
            LC::scalar en_BG;
        };


        // This is the information that will be extracted
        std::vector<InteractionPotential> interaction_data;
        //std::vector<RadialInteractionPotential> interaction_data;
        int thetaPoints = _widget.interactionThetaPoints;
        int phiPoints = _widget.interactionPhiPoints;

#define FULL_HELIKNOTON_INTERACTION 1

        // Spherical interaction
        if (radioInteractionType == 0) {


            for (int ti = 0; ti < thetaPoints; ti++) {
                bool degeneratePhi = false;
                for (int pi = 0; pi < phiPoints; pi++) {

                    LC::scalar theta, phi;
                    if (!_widget.interactionSymmetry) {
                        phi = (LC::scalar)pi / (phiPoints - 1) * 2. * M_PI;

                        if (thetaPoints > 1) // Don't waste time on degenerate polar points
                            theta = (LC::scalar)(ti + 1) / (thetaPoints + 1) * M_PI;
                        else
                            theta = M_PI * 0.5;
                    }
                    else {
                        phi = (LC::scalar)pi / (phiPoints - 1) * 2. * M_PI;
                        if (thetaPoints > 1) {
                            // Don't waste time on degenerate polar points
                            // If theta = 0, only add one point
                            theta = (LC::scalar)ti / (thetaPoints - 1) * 0.5 * M_PI;
                            if (ti == 0)
                                degeneratePhi = true;
                        }
                        else
                            theta = 0.5 * M_PI;
                    }

                    LC::scalar sph_r = 0.5 * _widget.separationDistance;

                    Eigen::Vector3d disp;
                    disp[0] = sph_r * sin(theta) * cos(phi);
                    disp[1] = sph_r * sin(theta) * sin(phi);
                    disp[2] = sph_r * cos(theta);

                    InteractionPotential interaction;
                    interaction.theta = theta;
                    interaction.phi = phi;
                    interaction.separation = 2.0 * disp.norm();


                    translations.emplace_back(disp);
                    interaction_data.push_back(interaction);

                    // Skip rest of phi iterations
                    if (degeneratePhi)
                        break;
                }
            }
        }
        else if (radioInteractionType == 1) { // Radial interaction

            bool degeneratePhi = false;
            for (int ri = 0; ri < r_points; ri++) {

                LC::scalar r;
                if (r_points == 1)
                    r = start_dist;
                else
                    r = start_dist + (LC::scalar)ri / (r_points - 1) * (end_dist - start_dist);


                LC::scalar rr = 0.5 * r;

                InteractionPotential interaction;
                interaction.theta = theta0;
                interaction.phi = phi0;
                interaction.separation = r;

                Eigen::Vector3d displacement;
                displacement[0] = rr * sin(theta0 * M_PI / 180.f) * cos(phi0 * M_PI / 180.f);
                displacement[1] = rr * sin(theta0 * M_PI / 180.f) * sin(phi0 * M_PI / 180.f);
                displacement[2] = rr * cos(theta0 * M_PI / 180.f);

                translations.emplace_back(displacement);
                interaction_data.push_back(interaction);

                // Skip rest of phi iterations
                if (degeneratePhi)
                    break;
            }

        }
        else if (radioInteractionType == 2) { // Volumetric interaction

            for (int ri = 0; ri < r_points; ri++) {

                LC::scalar r;
                if (r_points == 1)
                    r = start_dist;
                else
                    r = start_dist + (LC::scalar)ri / (r_points - 1) * (end_dist - start_dist);

                LC::scalar rr = 0.5 * r;

                for (int ti = 0; ti < thetaPoints; ti++) {
                    bool degeneratePhi = false;
                    for (int pi = 0; pi < phiPoints; pi++) {

                        LC::scalar theta, phi;
                        if (!_widget.interactionSymmetry) {
                            phi = (LC::scalar)pi / (phiPoints - 1) * 2. * M_PI;

                            if (thetaPoints > 1) // Don't waste time on degenerate polar points
                                theta = (LC::scalar)(ti + 1) / (thetaPoints + 1) * M_PI;
                            else
                                theta = M_PI * 0.5;

                        }
                        else {
                            phi = (LC::scalar)pi / (phiPoints - 1) * 2. * M_PI;
                            if (thetaPoints > 1) {
                                // Don't waste time on degenerate polar points
                                // If theta = 0, only add one point
                                theta = (LC::scalar)ti / (thetaPoints - 1) * 0.5 * M_PI;
                                if (ti == 0)
                                    degeneratePhi = true;
                            }
                            else
                                theta = 0.5 * M_PI;
                        }
                        // For testing purposes
#if FULL_HELIKNOTON_INTERACTION==0
                        theta = theta0;
                        phi = phi0;
#endif
                        Eigen::Vector3d disp;
                        disp[0] = rr * sin(theta) * cos(phi);
                        disp[1] = rr * sin(theta) * sin(phi);
                        disp[2] = rr * cos(theta);

                        InteractionPotential interaction;
                        interaction.theta = theta;
                        interaction.phi = phi;
                        interaction.separation = r;


                        translations.emplace_back(disp);
                        interaction_data.push_back(interaction);

                        // Skip rest of phi iterations
                        if (degeneratePhi)
                            break;
                    }
                }
            }
        }
        else {
            LC_ERROR("No translation points allocated");
        }

        // Iterate through translations, relaxing and computing energies

        // Get the zero point energy
        // Note that the total energy must be muiltiplied by dx, dy, or dz for the integration scaling to work out...
        // Compute the zero point energy

        LC::scalar T = 293; // 20 C in Kelvin
        // kb
        LC::scalar kb = 1.380649e-23;
        LC::scalar pitch = 1e-6; // um
        LC::scalar K = (data->k11.first + data->k22.first + data->k33.first) / 3. * 1e-12; // Avg elastic constant (Newtons)
        LC::scalar Kp = K * pitch;

        LC::scalar kbT_conversion = Kp / (kb * T);
        LC_INFO("Conversion factor = {0}", kbT_conversion);

        clear_heliknotons();
        solver->SetVoltage(EField * data->cell_dims[2], 50);
        LC::scalar E_bg = kbT_conversion * solver->TotalEnergy();
        LC::scalar E_H = 0.;
        LC::scalar z_prev = -1000000.;

        std::vector<uint32_t> full_region;
        uint32_t volNew = vNew[0] * vNew[1] * vNew[2];
        for (int i = 0; i < vNew[0] * vNew[1] * vNew[2]; i++) {
            full_region.push_back(i);
        }

        std::unique_ptr<float[]> preimage(new float[volNew]);
        LC::ExtendedMC::MarchingCubes mc;
        mc.clean_all();
        mc.set_resolution(vNew[0], vNew[1], vNew[2]);
        mc.init_all();
        mc.set_ext_data(preimage.get());

        // Run the test
#if FULL_HELIKNOTON_INTERACTION
        int count = 0;
        for (const auto& translation : translations) {

            if (abs(z_prev - translation.z()) > 1e-6) { // Recompute heliknoton energy for this layer
                clear_heliknotons();

                // Center the interaction in the middle of the volume
                embed_heliknoton(translation);

                // Set voltage and relax voltage
                solver->SetVoltage(EField * data->cell_dims[2], 50);
                solver->DomainRelax(_widget.interactionIterations, full_region, true, true);

                // Compute the energy for two separate heliknotons

                E_H = kbT_conversion * solver->TotalEnergy();
                z_prev = translation.z();
                LC_INFO("E_bg = {0}\tE_H = {1}", E_bg, E_H);

                // This is the entire calculation
                if (_widget.singleInteractionHeliknoton) {
                    // Compute and save the energy
                    interaction_data[count].energy = E_H - E_bg; // Single heliknoton confinement energy
                    interaction_data[count].en_HH = 0.;
                    interaction_data[count].en_H = E_H;
                    interaction_data[count].en_BG = E_bg;
                    interaction_data[count++].separation *= 0.5; // Need to multiply by half to get the distance from the center
                }
            }
            else if (_widget.singleInteractionHeliknoton) { // Same energy as before
                interaction_data[count].energy = E_H - E_bg; // Single heliknoton confinement energy
                interaction_data[count].en_HH = 0.;
                interaction_data[count].en_H = E_H;
                interaction_data[count].en_BG = E_bg;
                interaction_data[count++].separation *= 0.5; // Need to multiply by half to get the distance from the center
            }

            // Create full configuration
            if (!_widget.singleInteractionHeliknoton) {
                // Reset
                clear_heliknotons();
                embed_heliknoton(translation);
                embed_heliknoton(-translation);

                // Set voltage and relax voltage
                solver->SetVoltage(EField * data->cell_dims[2], 50);
                solver->DomainRelax(_widget.interactionIterations, full_region, true, true);

                LC::scalar E_HH = kbT_conversion * solver->TotalEnergy();


                // Compute and save the energy
                interaction_data[count].energy = E_HH - 2. * E_H + E_bg;// -zero_point_en;
                interaction_data[count].en_HH = E_HH;
                interaction_data[count].en_H = E_H;
                interaction_data[count].en_BG = E_bg;

                std::vector<Eigen::Vector3d> translation_pair;
                translation_pair.push_back(translation);
                translation_pair.push_back(-translation);


                LC_INFO("Energy: {2}, Iteration {0}/{1}", count + 1, translations.size(), interaction_data[count].energy);
                ++count;
            }
        }
        
        std::ofstream ofile(saveName.c_str(), std::ios::out | std::ios::binary);

        if (!ofile) {
            LC_WARN("Failed to save file <{0}>", saveName.c_str());
        }
        else {
            // Write file version
            char interactionType = 'e'; // Default: e -> error
            if (radioInteractionType == 0) {
                interactionType = 's';
            }
            else if (radioInteractionType == 1) {
                interactionType = 'r';
            }
            else if (radioInteractionType == 2) {
                interactionType = 'v';
            }
            else {
                LC_ERROR("Unknown interaction type [Invalid file]");
            }
            ofile.write((char*)&interactionType, sizeof(char));
            ofile.write((char*)&interaction_data[0], interaction_data.size() * sizeof(InteractionPotential));
        }

#else
            // See what the result is without doing the computation
        for (const auto& translation : translations) {
            embed_heliknoton(translation);
            embed_heliknoton(-translation);
        }
#endif

        initVisuals();
        _widget.updateImage = true;
    }

    ImGui::End();
}


void Sandbox::findVortexKnotComponents(bool saveObj, const std::string& obj_name) {
    using namespace VKH;

    Dataset* data = static_cast<Dataset*>(_solver->GetDataPtr());
    std::array<float, 3> cellf = {
        static_cast<float>(data->cell_dims[0]),
        static_cast<float>(data->cell_dims[1]),
        static_cast<float>(data->cell_dims[2])
    };

    GridMetrics base_grid(data->voxels, cellf, 1);
    GridMetrics interp_grid(data->voxels, cellf, _widget.knot_interaction_handle.upsampling_multiplier);

    for (auto& knot : _vortexKnot) {
        knot.draw = false;
    }
    _vortexKnot.clear();
    _vortexManipulator = std::make_unique<LC::Drawable::Object3D>();
    _vortexManipulator->setParent(_manipulator.get());

    auto smoothMCMesh = [&](LC::ExtendedMC::MarchingCubes& mc_inst, float lambda, float mu, int iterations, float* nn_ptr, const std::array<float, 4>& color) -> SmoothedMesh {
        if (mc_inst.nverts() == 0 || mc_inst.ntrigs() == 0) return {};

        int nverts = mc_inst.nverts();
        auto iso_verts = std::make_unique<LC::Math::IsoVertex[]>(nverts);
        LC::ExtendedMC::Vertex* raw_verts = mc_inst.vertices();

        for (int i = 0; i < nverts; ++i) {
            iso_verts[i] = toIsoVertex(raw_verts[i]);
        }

        uint* indices = reinterpret_cast<uint*>(mc_inst.triangles());
        auto result = TaubinSmoothingGeneral(
            lambda, mu, iterations,
            iso_verts.get(), nverts, indices, mc_inst.ntrigs() * 3,
            nn_ptr, data->voxels, color, 0, 0, 0, 0
        );

        return { std::move(std::get<0>(result)), std::move(std::get<1>(result)) };
    };

    auto addKnotToScene = [&](SmoothedMesh& mesh, const Color4& color, bool is_ribbon) {
        _vortexKnot.emplace_back();
        VortexKnot& knot = _vortexKnot.back();
        knot.knotColor = color;
        knot.ribbon = is_ribbon;

        knot.surface.Init(
            reinterpret_cast<LC::Surface::Vertex*>(mesh.vertices.data()),
            mesh.vertices.size(),
            &mesh.triangles[0][0],
            mesh.triangles.size() * 3,
            _widget.preimage_translate
        );
        knot.mesh = knot.surface.Mesh();

        new LC::Drawable::TransparentNormalDrawable{
            *_vortexManipulator, _phongShader, *(knot.mesh),
            knot.draw, knot.cullFaces, _transparentNormalDrawables
        };
    };

    std::unique_ptr<float[]> chi_field;
    auto field_nn = std::make_unique<float[]>(3 * base_grid.vol);
    auto valid_fieldf = std::make_unique<float[]>(base_grid.vol);
    auto field_S = std::make_unique<float[]>(base_grid.vol);

    for (uint id = 0; id < 3 * base_grid.vol; ++id) {
        field_nn[id] = static_cast<float>(data->directors[id]);
    }

    LC::Math::MeshProcessor processor( base_grid.dims[0], base_grid.dims[1], base_grid.dims[2] );
    processor.SetSpacing(base_grid.dr[0], base_grid.dr[1], base_grid.dr[2]);
    processor.SetOrigin( -0.5f * cellf[0], -0.5f * cellf[1], -0.5f * cellf[2] );

    float vortexTubeRadius = _vortex_line.tubeRadius * 0.35f * (std::min)({ base_grid.dr[0], base_grid.dr[1], base_grid.dr[2] });
    auto umbilicMesh = processor.ComputeUmbilicTube(field_nn.get(), vortexTubeRadius, 12, 0.2f);

    _vortex_line.numNodes = LC::Math::ChiralityField(
        field_nn.get(), chi_field, data->voxels, cellf,
        _vortex_line.valid_field, true, _widget.knot_interaction_handle.tolerance, false
    );

    computeOrderParameterAndFilter(
        base_grid, field_nn.get(), chi_field.get(),
        _widget.knot_interaction_handle.isoValue,
        _vortex_line.valid_field.get(), valid_fieldf.get(), field_S.get()
    );

    _vortex_line.numNodes = 0;
    for (uint ii = 0; ii < base_grid.vol; ++ii) {
        if (valid_fieldf[ii] < 0.f) {
            _vortex_line.numNodes++;
        }
    }


    std::array<float, 4> col_white = { 1.f, 1.f, 1.f, 1.f };
    auto sample_grid = std::make_unique<float[]>(interp_grid.vol);
    std::vector<uint8_t> upsampled_mask(interp_grid.vol, 0);
    std::vector<uint> vortex_indices_upsampled;

    LC::ExtendedMC::MarchingCubes mc2;
    mc2.clean_all();
    mc2.set_resolution(interp_grid.dims[0], interp_grid.dims[1], interp_grid.dims[2]);
    mc2.set_ext_data(sample_grid.get());
    mc2.init_all();

    if (!umbilicMesh.empty()) {
        // Format the mesh
        std::vector<MeshLib::PNCVertex<float>> umbilic_vertices;
        std::vector<MeshLib::Triangle> umbilic_triangles;
        size_t nVert = umbilicMesh.numVertices();
        size_t nTrig = umbilicMesh.numTriangles();
        umbilic_vertices.resize(nVert);
        umbilic_triangles.resize(nTrig);

        float normSign = 1.0f;
        auto color = Magnum::Color4{ 0.5f, 0.5f, 0.5f, 0.5f };

        for (size_t i = 0; i < nVert; ++i) {
            const auto& src = umbilicMesh.vertices[i];
            auto& dst = umbilic_vertices[i];

            dst.position[0] = src.position[0];
            dst.position[1] = src.position[1];
            dst.position[2] = src.position[2];

            dst.normal[0] = normSign * src.normal[0];
            dst.normal[1] = normSign * src.normal[1];
            dst.normal[2] = normSign * src.normal[2];

            dst.color[0] = color[0];
            dst.color[1] = color[1];
            dst.color[2] = color[2];
            dst.color[3] = color[3];
        }

        for (size_t t = 0; t < nTrig; ++t) {
            umbilic_triangles[t][0] = umbilicMesh.indices[3 * t + 0];
            umbilic_triangles[t][1] = umbilicMesh.indices[3 * t + 2];
            umbilic_triangles[t][2] = umbilicMesh.indices[3 * t + 1];
        }
        _vortex_line.components = extractConnectedComponents(SmoothedMesh{umbilic_vertices, umbilic_triangles});

        std::array<LC::scalar, 3> cell = {
            data->cell_dims[0] / (data->voxels[0] - 1),
            data->cell_dims[1] / (data->voxels[1] - 1),
            data->cell_dims[2] / (data->voxels[2] - 1)
        };

        int col_i = 0;
        float d_deg = _vortex_line.components.empty() ? 0.f : 360.f / static_cast<float>(_vortex_line.components.size());

        for (const auto& component : _vortex_line.components) {
            Color4 color = Color3::fromHsv({ Deg(d_deg * col_i++), 1.f, 1.f });
            std::array<float, 4> col_arr = { color[0], color[1], color[2], color[3] };

            // Get the subset of vertices and relabel triangles
            std::vector<MeshLib::PNCVertex<float>> comp_vertices(component.size());
            std::vector<MeshLib::Triangle> comp_triangles;
            
            comp_triangles.reserve(2*component.size());

            std::map<size_t, size_t> map_glob_to_local;

            for (size_t i = 0; i < component.size(); ++i) {
                size_t c_i = component[i];
                comp_vertices[i] = umbilic_vertices[c_i];
                map_glob_to_local.emplace(c_i, i);
            }

            for (const auto& tri : umbilic_triangles) {
                // Look for the first index to see if the triangle will be in the list
                auto it = map_glob_to_local.find(tri[0]);
                if (it != map_glob_to_local.end()) {
                    MeshLib::Triangle local_triangle = { map_glob_to_local[tri[0]], map_glob_to_local[tri[1]], map_glob_to_local[tri[2]] };
                    comp_triangles.emplace_back(local_triangle);
                }
            }
           
            colorVerticesByDirector(comp_vertices, base_grid, interp_grid, field_nn.get(), upsampled_mask, vortex_indices_upsampled);
            addKnotToScene(SmoothedMesh{ comp_vertices, comp_triangles }, color, false);
        }

        LC_INFO("Discovered {0} component(s)", _vortex_line.components.size());
        for (const auto& c : _vortex_line.components) {
            LC_INFO("- Component size = {0}", c.size());
        }
    }

    auto field_Q = std::make_unique<float[]>(5 * base_grid.vol);
    computeChiQTensor(base_grid, chi_field.get(), field_Q.get());

    std::vector<uint> fat_vortex_list = buildFatVortexRegion(
        interp_grid, vortex_indices_upsampled, upsampled_mask,
        _widget.knot_interaction_handle.sb_saturation
    );

    auto field_Ssb = std::make_unique<float[]>(interp_grid.vol);
    SplayBendStats sb_stats = computeInterpolatedSplayBend(
        base_grid, interp_grid, field_Q.get(), fat_vortex_list, field_Ssb.get()
    );

    auto buildSplayBendRibbon = [&](float threshold_ratio, float extremum_sb, float lambda, const std::array<float, 4>& col_arr, const Color4& col4, const char* label) {
        populateSplayBendSampleGrid(interp_grid, fat_vortex_list, field_Ssb.get(), threshold_ratio, extremum_sb, sample_grid.get());
        mc2.reset_mesh();
        mc2.run();
        rescaleAndNormalizeMC(mc2, interp_grid);

        SmoothedMesh ribbon_mesh = smoothMCMesh(mc2, lambda, -0.34f, 80, field_nn.get(), col_arr);
        mc2.reset_mesh();

        if (!ribbon_mesh.empty()) {
            addKnotToScene(ribbon_mesh, col4, true);
        }
        else {
            LC_INFO("Not enough points for {0}", label);
        }
    };

    buildSplayBendRibbon(sb_stats.tc_pos, sb_stats.max_sb, 0.35f, { 0.f, 0.f, 1.f, 1.f }, Color4{ 0.f, 0.f, 1.f, 1.f }, "positive SB");
    buildSplayBendRibbon(sb_stats.tc_neg, sb_stats.min_sb, 0.33f, { 1.f, 1.f, 0.f, 1.f }, Color4{ 1.f, 1.f, 0.f, 1.f }, "negative SB");
}

std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>>
Sandbox::findVortexKnot(
    bool saveObj,
    const std::string& obj_name,
    std::unique_ptr<float[]>& chi_field,
    std::unique_ptr<float[]>& field_nn,
    std::unique_ptr<float[]>& valid_fieldf,
    std::unique_ptr<float[]>& field_S,
    std::unique_ptr<float[]>& sample_grid,
    std::unique_ptr<float[]>& field_Q,
    std::unique_ptr<float[]>& field_Ssb,
    LC::ExtendedMC::MarchingCubes& mc,
    LC::ExtendedMC::MarchingCubes& mc2,
    std::array<int, 3>& vox_interp)
{
    using namespace VKH;
    Dataset* data = static_cast<Dataset*>(_solver->GetDataPtr());
    std::array<float, 3> cellf = {
        static_cast<float>(data->cell_dims[0]),
        static_cast<float>(data->cell_dims[1]),
        static_cast<float>(data->cell_dims[2])
    };

    GridMetrics base_grid(data->voxels, cellf, 1);
    GridMetrics interp_grid(vox_interp, cellf, 1);

    for (auto& knot : _vortexKnot) {
        knot.draw = false;
    }
    _vortexKnot.clear();

    for (uint id = 0; id < 3 * base_grid.vol; ++id) {
        field_nn[id] = static_cast<float>(data->directors[id]);
    }

    LC::Math::ChiralityField(
        field_nn.get(), chi_field, data->voxels, cellf,
        _vortex_line.valid_field, false, _widget.knot_interaction_handle.tolerance, false
    );

    computeOrderParameterAndFilter(
        base_grid, field_nn.get(), chi_field.get(),
        _widget.knot_interaction_handle.isoValue,
        _vortex_line.valid_field.get(), valid_fieldf.get(), field_S.get(),
        true
    );

    _vortex_line.numNodes = 0;
    for (uint ii = 0; ii < base_grid.vol; ++ii) {
        if (valid_fieldf[ii] < 0.f) {
            _vortex_line.numNodes++;
        }
    }

    auto smoothAndExportMC = [&](
        LC::ExtendedMC::MarchingCubes& mc_inst,
        float lambda, float mu, int iterations,
        float* nn_ptr, const std::array<float, 4>& color,
        int smooth_mode, bool do_save, const std::string& filename,
        std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>* out_raw_mesh = nullptr) -> SmoothedMesh
        {
            if (mc_inst.nverts() == 0 || mc_inst.ntrigs() == 0) return {};

            int nverts = mc_inst.nverts();
            int ntrigs = mc_inst.ntrigs();
            auto iso_verts = std::make_unique<LC::Math::IsoVertex[]>(nverts);
            LC::ExtendedMC::Vertex* raw_verts = mc_inst.vertices();
            uint* indices = reinterpret_cast<uint*>(mc_inst.triangles());

            for (int i = 0; i < nverts; ++i) {
                iso_verts[i] = toIsoVertex(raw_verts[i]);
            }

            if (out_raw_mesh) {
                auto& out_v = std::get<0>(*out_raw_mesh);
                auto& out_i = std::get<1>(*out_raw_mesh);
                out_v.assign(iso_verts.get(), iso_verts.get() + nverts);
                out_i.assign(indices, indices + ntrigs * 3);
            }

            if (do_save || !filename.empty()) {
                auto result = TaubinSmoothingGeneral(
                    lambda, mu, iterations,
                    iso_verts.get(), nverts, indices, ntrigs * 3,
                    nn_ptr, data->voxels, color, false, smooth_mode, 0, do_save, filename
                );
                return { std::move(std::get<0>(result)), std::move(std::get<1>(result)) };
            }
            else {
                auto result = TaubinSmoothingGeneral(
                    lambda, mu, iterations,
                    iso_verts.get(), nverts, indices, ntrigs * 3,
                    nn_ptr, data->voxels, color, 0, 0, 0, 0
                );
                return { std::move(std::get<0>(result)), std::move(std::get<1>(result)) };
            }
    };

    mc.set_ext_data(valid_fieldf.get());
    mc2.set_ext_data(sample_grid.get());
    mc.reset_mesh();
    mc2.reset_mesh();
    mc.run();

    std::array<float, 4> col_white = { 1.f, 1.f, 1.f, 1.f };
    SmoothedMesh initial_mesh = smoothAndExportMC(mc, 0.5f, -0.34f, 80, nullptr, col_white, 0, false, "");
    mc.reset_mesh();
    rescaleAndNormalizeMesh(initial_mesh, base_grid);

    if (!initial_mesh.empty()) {
        std::vector<std::vector<uint>> components = extractConnectedComponents(initial_mesh);

        if (_widget.knot_interaction_handle.saveStringKnots) {
            LC::scalar dist = 2.0 * _widget.knot_interaction_handle.point_density *
                std::sqrt(base_grid.dr[0] * base_grid.dr[0] + base_grid.dr[1] * base_grid.dr[1] + base_grid.dr[2] * base_grid.dr[2]);

            auto exported_data = exportVortexKnot(initial_mesh.vertices, initial_mesh.triangles, dist, obj_name + "_vortex.bin");
            const auto& pdata = std::get<0>(exported_data);
            std::vector<Eigen::Vector2d> odata(pdata.size());

            for (size_t p = 0; p < pdata.size(); ++p) {
                int i = static_cast<int>((pdata[p].x() / cellf[0] + 0.5) * (base_grid.dims[0] - 1));
                int j = static_cast<int>((pdata[p].y() / cellf[1] + 0.5) * (base_grid.dims[1] - 1));
                int k = static_cast<int>((pdata[p].z() / cellf[2] + 0.5) * (base_grid.dims[2] - 1));

                if (base_grid.inBounds(i, j, k)) {
                    uint id = base_grid.to1D(i, j, k);
                    odata[p][0] = std::acos(field_nn[id + 2 * base_grid.vol]);
                    odata[p][1] = std::atan2(field_nn[id + base_grid.vol], field_nn[id]);
                }
            }

            std::ofstream ofile(obj_name + "_sphere.bin", std::ios::out | std::ios::binary);
            if (ofile.is_open()) {
                int32_t osz = static_cast<int32_t>(odata.size());
                ofile.write(reinterpret_cast<const char*>(&osz), sizeof(int32_t));
                ofile.write(reinterpret_cast<const char*>(odata.data()), osz * sizeof(Eigen::Vector2d));
            }
        }

        int radius = _widget.knot_interaction_handle.saturation;
        int comp_idx = 0;

        for (const auto& component : components) {
            stampComponentToGrid(component, initial_mesh.vertices, interp_grid, radius, sample_grid.get());

            mc2.reset_mesh();
            mc2.run();
            rescaleAndNormalizeMC(mc2, interp_grid);

            std::string ply_name = obj_name + "_vortex" + std::to_string(++comp_idx) + "_of_" + std::to_string(components.size()) + ".ply";
            smoothAndExportMC(mc2, 0.5f, -0.34f, 120, field_nn.get(), col_white, 1, true, ply_name);
        }
    }

#define SAVE_SPLAY_BEND 1
#if SAVE_SPLAY_BEND
    computeChiQTensor(base_grid, chi_field.get(), field_Q.get());

    float stdev = computeBaseGridSplayBend(
        base_grid, valid_fieldf.get(), field_Q.get(),
        _widget.knot_interaction_handle.sb_saturation, field_Ssb.get()
    );

    auto exportSplayBendBand = [&](float lower_bound, float upper_bound, float lambda, const std::array<float, 4>& color, const std::string& suffix) {
        std::fill_n(valid_fieldf.get(), base_grid.vol, 0.f);

        for (int k = 1; k < base_grid.dims[2] - 1; ++k) {
            for (int j = 1; j < base_grid.dims[1] - 1; ++j) {
                for (int i = 1; i < base_grid.dims[0] - 1; ++i) {
                    uint id = base_grid.to1D(i, j, k);
                    if (lower_bound < field_Ssb[id] && field_Ssb[id] < upper_bound) {
                        valid_fieldf[id] = -1.f;
                    }
                }
            }
        }

        mc.reset_mesh();
        mc.run();
        rescaleAndNormalizeMC(mc, base_grid);

        smoothAndExportMC(mc, lambda, -0.34f, 80, field_nn.get(), color, 0, saveObj, obj_name + suffix);
        mc.reset_mesh();
        };

    exportSplayBendBand(
        _widget.knot_interaction_handle.sb_pos_iso_ratio_lower * stdev,
        _widget.knot_interaction_handle.sb_pos_iso_ratio_upper * stdev,
        0.35f, { 0.f, 0.f, 1.f, 1.f }, "_SB_pos.ply"
    );

    exportSplayBendBand(
        -_widget.knot_interaction_handle.sb_neg_iso_ratio_upper * stdev,
        -_widget.knot_interaction_handle.sb_neg_iso_ratio_lower * stdev,
        0.33f, { 1.f, 1.f, 0.f, 1.f }, "_SB_neg.ply"
    );
#endif

    int preim_ct = 0;
    int nPreimages = static_cast<int>(_widget.knot_interaction_handle.preimages.size());

    std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>> heliknoton_list;
    heliknoton_list.reserve(nPreimages);

    for (const auto& preim : _widget.knot_interaction_handle.preimages) {
        float iso = _widget.knot_interaction_handle.preimage_isovalue;
        float th = preim.theta / 180.f * static_cast<float>(M_PI);
        float ph = preim.phi / 180.f * static_cast<float>(M_PI);

        Eigen::Vector3f p0(std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), std::cos(th));

        populatePreimageSampleGrid(base_grid, interp_grid, field_nn.get(), p0, sample_grid.get());

        mc2.reset_mesh();
        mc2.run(iso);
        rescaleAndNormalizeMC(mc2, interp_grid);

        if (mc2.nverts() > 0 && mc2.ntrigs() > 0) {
            auto col = LC::Imaging::Colors::RungeSphere(th, 2.f * ph);
            std::array<float, 4> col4 = { col[0], col[1], col[2], 1.f };

            std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>> p_list;
            std::string fullName = obj_name + "_preim_" + std::to_string(++preim_ct) + "_of_" + std::to_string(nPreimages) + ".ply";

            smoothAndExportMC(mc2, 0.33f, -0.34f, 80, field_nn.get(), col4, 0, saveObj, fullName, &p_list);
            mc2.reset_mesh();

            heliknoton_list.emplace_back(std::move(p_list));
        }
    }

    return heliknoton_list;
}

std::tuple <std::vector<MeshLib::PNCVertex<float>>, std::vector<MeshLib::Triangle>>
Sandbox::TaubinSmoothingGeneral(
    float lambda,
    float mu,
    int smoothingIterations,
    LC::Math::IsoVertex* verts,
    unsigned int nVert,
    unsigned int* indices,
    unsigned int nInd,
    const float* field_nn,
    const std::array<int, 3>& Vox,
    const std::array<float, 4>& color,
    bool invert_normals,
    int color_style,
    unsigned int nSubdivisions,
    bool saveObj,
    const std::string& obj_name) {

    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    unsigned int size = Vox[0] * Vox[1] * Vox[2];

    // Subdivide the mesh
    LC::Algorithm::Mesh<float> mesh;
    std::vector<MeshLib::PNCVertex<float>> vertices;
    std::vector<MeshLib::Triangle> triangles;

    mesh.read_pnc_mesh((MeshLib::PNCVertex<float>*)verts, nVert, (MeshLib::Triangle*)indices, nInd / 3);
    LC::Algorithm::loop_subdivision(&mesh, nSubdivisions);
    // Export mesh data to vertices and indices (Don't compute normals yet)
    mesh.write_obj(vertices, triangles, false);

    // Start Taubin smoothing (https://graphics.stanford.edu/courses/cs468-01-fall/Papers/taubin-smoothing.pdf)
    {
        // Use graph to find edges
        Graph graph(vertices.size());
        std::vector<MeshLib::PNCVertex<float>> vertices_prime;
        vertices_prime.resize(vertices.size());

        for (const auto& tri : triangles) {
            // Each face contributes 3 edges for a triangle
            for (int d = 0; d < 3; d++) {
                uint dp1 = (d + 1) % 3;
                graph.addEdge(tri[d], tri[dp1]);
            }
        }

        // Weight
        float w_ij;
        // Smoothing parameters (0 < lambda < -mu)
        std::list<uint>* adjacencyList = graph.adjacencyList();

        for (int s = 0; s < smoothingIterations; s++) {

            // Positive Gaussian smoothing step
            for (int i = 0; i < vertices.size(); i++) {
                Eigen::Vector3f v_i(vertices[i].position[0], vertices[i].position[1], vertices[i].position[2]);

                // Sum up to second neighbors average displacement
                Eigen::Vector3f v_avg(0., 0., 0.);
                w_ij = 1. / adjacencyList[i].size();
                for (auto j : adjacencyList[i]) {
                    Eigen::Vector3f v_j(vertices[j].position[0], vertices[j].position[1], vertices[j].position[2]);
                    // Nearest neighbors
                    //v_avg = v_avg + w_ij * (v_j - v_i);
                    // Second nearest neighbors
                    float w_ijk = w_ij / adjacencyList[j].size();
                    for (auto k : adjacencyList[j]) {
                        Eigen::Vector3f v_k(vertices[k].position[0], vertices[k].position[1], vertices[k].position[2]);
                        v_avg = v_avg + w_ijk * (v_k - v_i);
                    }

                }
                Eigen::Vector3f vp = v_i + lambda * v_avg;
                vertices_prime[i] = { vp(0), vp(1), vp(2) };
            }

            // Update vertices
            for (int i = 0; i < vertices.size(); i++)
                vertices[i].position = vertices_prime[i].position;

            // Negative Gaussian smoothing step
            for (int i = 0; i < vertices.size(); i++) {
                Eigen::Vector3f v_i(vertices[i].position[0], vertices[i].position[1], vertices[i].position[2]);

                Eigen::Vector3f v_avg(0., 0., 0.);
                // Sum over neighbors average displacement
                w_ij = 1. / adjacencyList[i].size();
                for (auto j : adjacencyList[i]) {
                    Eigen::Vector3f v_j(vertices[j].position[0], vertices[j].position[1], vertices[j].position[2]);
                    //v_avg = v_avg + w_ij * (v_j - v_i);
                    // Second nearest neighbors
                    float w_ijk = w_ij / adjacencyList[j].size();
                    for (auto k : adjacencyList[j]) {
                        Eigen::Vector3f v_k(vertices[k].position[0], vertices[k].position[1], vertices[k].position[2]);
                        v_avg = v_avg + w_ijk * (v_k - v_i);
                    }
                }
                Eigen::Vector3f vp = v_i + mu * v_avg;
                vertices_prime[i] = { vp(0), vp(1), vp(2) };
            }

            // Update vertices
            for (int i = 0; i < vertices.size(); i++)
                vertices[i].position = vertices_prime[i].position;

        }
    }

    // Compute normals
    MeshLib::compute_normals(vertices, triangles, invert_normals);

    std::array<float, 3> cell = { (float)data->cell_dims[0], (float)data->cell_dims[1], (float)data->cell_dims[2] };

    // Set the color
    for (auto& v : vertices) {
        if (color_style == 1) { // Pontryagin
            // Color depends on orientation
            int ii = (v.position[0] / data->cell_dims[0] + 0.5) * (Vox[0] - 1);
            int jj = (v.position[1] / data->cell_dims[1] + 0.5) * (Vox[1] - 1);
            int kk = (v.position[2] / data->cell_dims[2] + 0.5) * (Vox[2] - 1);

            uint id = ii + jj * Vox[0] + kk * Vox[0] * Vox[1];
            float theta = acos(field_nn[id + 2 * size]);
            float phi = atan2(field_nn[id + size], field_nn[id]);
            auto col = LC::Imaging::Colors::RungeSphere(theta, phi);
            for (int d = 0; d < 3; d++)
                v.color[d] = col[d];

            v.color[3] = color[3];
        }
        else if (color_style == 2) { // Pontryagin (for S2/Z2)
            // Color depends on orientation
            int ii = (v.position[0] / data->cell_dims[0] + 0.5) * (Vox[0] - 1);
            int jj = (v.position[1] / data->cell_dims[1] + 0.5) * (Vox[1] - 1);
            int kk = (v.position[2] / data->cell_dims[2] + 0.5) * (Vox[2] - 1);

            uint id = ii + jj * Vox[0] + kk * Vox[0] * Vox[1];
            float theta = acos(field_nn[id + 2 * size]);
            float phi = atan2(field_nn[id + size], field_nn[id]);
            auto col = LC::Imaging::Colors::RungeSphere(theta, 2.f * phi);
            for (int d = 0; d < 3; d++)
                v.color[d] = col[d];

            v.color[3] = color[3];
        }
        else if (color_style == 3) { // Helicity
            // Color depends on orientation
            int ii = (v.position[0] / data->cell_dims[0] + 0.5) * (Vox[0] - 1);
            int jj = (v.position[1] / data->cell_dims[1] + 0.5) * (Vox[1] - 1);
            int kk = (v.position[2] / data->cell_dims[2] + 0.5) * (Vox[2] - 1);

            uint id = ii + jj * Vox[0] + kk * Vox[0] * Vox[1];
            auto chi = LC::Math::HandednessTensor(ii, jj, kk, field_nn, Vox, cell);
            float result = -0.5f * chi.trace() / (2.f * M_PI);
            //Eigen::Matrix3d chig = Eigen::Matrix3d::Zero();
            //chig(2,2) = -2. * M_PI;
            //auto chidiff = (chi - chig)/(2.*M_PI);
            //float result = (chidiff * chidiff).trace();

            Magnum::Color3 col;
            if (result > 0.f) col = col = Magnum::Color3(result, 0.5f, 0.5f);
            else col = Magnum::Color3(0.5f, 0.5f, -result);


            for (int d = 0; d < 3; d++) {
                v.color[d] = col[d];
                if (v.color[d] >= 1.f) { // Reached color saturation
                    v.color[d] = 1.;
                }
            }

            v.color[3] = color[3];
        }
        else v.color = color;
    }

    if (saveObj) {
        LC::Algorithm::Mesh<float> mesh_t;
        mesh_t.read_pnc_mesh(vertices, triangles);
        mesh_t.write_ply(obj_name.c_str());
    }

    return { vertices, triangles };
}

std::tuple <std::vector<MeshLib::PNCVertex<float>>, std::vector<MeshLib::Triangle>>
    Sandbox::TaubinSmoothing(LC::Math::IsoVertex* verts,
    unsigned int nVert,
    unsigned int* indices,
    unsigned int nInd,
    const float* field_nn,
    const std::array<int, 3> &Vox,
    const std::array<float, 4>& color,
    bool invert_normals,
    int color_style,
    unsigned int nSubdivisions,
    bool saveObj,
    const std::string &obj_name) {

    return TaubinSmoothingGeneral(_widget.smoothingLambda, _widget.smoothingMu, _widget.smoothingIterations,
        verts, nVert, indices, nInd, field_nn, Vox, color, invert_normals, color_style, nSubdivisions, saveObj, obj_name);

}

void Sandbox::generateIsosurface() {
    Dataset* data = static_cast<Dataset*>(_solver->GetDataPtr());
    std::array<float, 3> cellf = {
        static_cast<float>(data->cell_dims[0]),
        static_cast<float>(data->cell_dims[1]),
        static_cast<float>(data->cell_dims[2])
    };
    std::array<LC::scalar, 3> cell = {
        data->cell_dims[0] / (data->voxels[0] - 1),
        data->cell_dims[1] / (data->voxels[1] - 1),
        data->cell_dims[2] / (data->voxels[2] - 1)
    };

    std::array<int, 3> vNew;
    for (int d = 0; d < 3; ++d) {
        vNew[d] = _widget.shrink_interval_end[d] - _widget.shrink_interval_begin[d] + 1;
        _widget.preimage_translate[d] = 0.5f * (_widget.shrink_interval_end[d] - data->voxels[d] + _widget.shrink_interval_begin[d] - 1.0f) * cell[d];
    }

    VKH::GridMetrics base_grid(data->voxels, cellf, 1);
    VKH::GridMetrics sub_grid(vNew, cellf, 1);

    auto field_nn = std::make_unique<float[]>(3 * sub_grid.vol);
    FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
    Eigen::TensorMap<Eigen::Tensor<float, 4>> nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);

    for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; ++k) {
        int kp = k - _widget.shrink_interval_begin[2] + 1;
        for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; ++j) {
            int jp = j - _widget.shrink_interval_begin[1] + 1;
            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; ++i) {
                int ip = i - _widget.shrink_interval_begin[0] + 1;
                for (int d = 0; d < 3; ++d) {
                    nn_new(ip, jp, kp, d) = nn(i, j, k, d);
                }
            }
        }
    }

    _preimageManipulator = std::make_unique<LC::Drawable::Object3D>();
    _preimageManipulator->setParent(_manipulator.get());

    auto registerSurfaceDrawable = [&](auto& target_obj, LC::Surface::Vertex* verts, size_t nVert, uint* indices, size_t nInd, bool& cull_flag) {
        target_obj.surface.Init(verts, nVert, indices, nInd, _widget.preimage_translate);
        target_obj.mesh = target_obj.surface.Mesh();
        new LC::Drawable::TransparentNormalDrawable{
            *_preimageManipulator, _phongShader, *(target_obj.mesh),
            target_obj.draw, cull_flag, _transparentNormalDrawables
        };
    };

    auto field = std::make_unique<float[]>(sub_grid.vol);
    LC::Math::MeshProcessor meshProc(vNew[0], vNew[1], vNew[2]);

    float dx = static_cast<float>(cell[0]);
    float dy = static_cast<float>(cell[1]);
    float dz = static_cast<float>(cell[2]);
    meshProc.SetSpacing(dx, dy, dz);
    meshProc.SetOrigin(
        -0.5f * (vNew[0] - 1) * dx,
        -0.5f * (vNew[1] - 1) * dy,
        -0.5f * (vNew[2] - 1) * dz
    );

    for (S2Fiber& pimage : _nematicPreimages) {
        float theta = pimage.theta / 180.f * static_cast<float>(M_PI);
        float phi = pimage.phi / 180.f * static_cast<float>(M_PI);
        Eigen::Vector3f N(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));

        auto col = LC::Imaging::Colors::RungeSphere(theta, phi);
        std::array<float, 4> color = { col[0], col[1], col[2], pimage.alpha };

        LC::Math::MeshProcessor::ExtractedMesh extracted;

        // Primary Method: Orthonormal basis (b1, b2) intersection for vector preimages
        if (!pimage.pontryagin) {
            std::array<float, 3> n0 = { N.x(), N.y(), N.z() };
            // Scale tube radius to the physical voxel spacing and isoLevel
            float tubeRadius = std::sqrt(pimage.isoLevel) * dx;
            extracted = meshProc.ComputeBasisPreimageTube(field_nn.get(), n0, tubeRadius, 12, 0.0f);
        }

        bool fallback_extraction = false;

        // Fallback Method: VTK Marching Cubes Isosurface (for Pontryagin or if basis intersection is empty)
        if (extracted.empty()) {
            for (uint i = 0; i < sub_grid.vol; ++i) {
                float d_z = field_nn[i + 2 * sub_grid.vol] - N.z();
                if (pimage.pontryagin) {
                    field[i] = d_z * d_z;
                }
                else {
                    float d_x = field_nn[i] - N.x();
                    float d_y = field_nn[i + sub_grid.vol] - N.y();
                    field[i] = d_x * d_x + d_y * d_y + d_z * d_z;
                }
            }

            VKH::applyDomainInversion(field.get(), sub_grid.vol, pimage.isoLevel, pimage.domain_inversion);
            extracted = meshProc.ComputeIsosurface(field.get(), pimage.isoLevel);
            fallback_extraction = true;
        }

        if (!extracted.empty()) {
            std::vector<MeshLib::PNCVertex<float>> vertices;
            std::vector<MeshLib::Triangle> triangles;

            if (fallback_extraction) {
                uint nVert = static_cast<uint>(extracted.numVertices());
                uint nInd = static_cast<uint>(extracted.indices.size());
                auto smoothing_result = TaubinSmoothing(
                    extracted.vertices.data(), nVert, extracted.indices.data(), nInd, field_nn.get(),
                    vNew, color, pimage.normal_inversion, pimage.pontryagin, pimage.subdivisions
                );
                vertices = std::move(std::get<0>(smoothing_result));
                triangles = std::move(std::get<1>(smoothing_result));
                LC_INFO("Successfully generated fallback surface (verts = {0}, triangles = {1})", vertices.size(), triangles.size());
            }
            else {
                size_t nVert = extracted.numVertices();
                size_t nTrig = extracted.numTriangles();
                vertices.resize(nVert);
                triangles.resize(nTrig);

                float normSign = pimage.normal_inversion ? -1.0f : 1.0f;

                for (size_t i = 0; i < nVert; ++i) {
                    const auto& src = extracted.vertices[i];
                    auto& dst = vertices[i];

                    dst.position[0] = src.position[0];
                    dst.position[1] = src.position[1];
                    dst.position[2] = src.position[2];

                    dst.normal[0] = normSign * src.normal[0];
                    dst.normal[1] = normSign * src.normal[1];
                    dst.normal[2] = normSign * src.normal[2];

                    dst.color[0] = color[0];
                    dst.color[1] = color[1];
                    dst.color[2] = color[2];
                    dst.color[3] = color[3];
                }

                for (size_t t = 0; t < nTrig; ++t) {
                    triangles[t][0] = extracted.indices[3 * t + 0];
                    if (pimage.normal_inversion) {
                        triangles[t][1] = extracted.indices[3 * t + 2];
                        triangles[t][2] = extracted.indices[3 * t + 1];
                    }
                    else {
                        triangles[t][1] = extracted.indices[3 * t + 1];
                        triangles[t][2] = extracted.indices[3 * t + 2];
                    }
                }

                LC_INFO("Successfully generated basis tube surface (verts = {0}, triangles = {1})", vertices.size(), triangles.size());
            }

            registerSurfaceDrawable(
                pimage, reinterpret_cast<LC::Surface::Vertex*>(vertices.data()),
                vertices.size(), &triangles[0][0], triangles.size() * 3, pimage.cullFaces
            );
        }
    }

    std::unique_ptr<float[]> chi_field;
    std::unique_ptr<short[]> valid_field;
    if (_widget.chiColorScheme != 2 && _vortexShell) {
        LC::Math::ChiralityField(field_nn.get(), chi_field, vNew, cellf, valid_field, true, _widget.isoLevel);
    }

    if (_widget.generateKnots) {
        for (auto& knot : _processedVortexKnot) knot.draw = false;
        _processedVortexKnot.clear();

        _processedVortexManipulator = std::make_unique<LC::Drawable::Object3D>();
        _processedVortexManipulator->setParent(_manipulator.get());

        if (_vortex_line.numNodes && !_vortex_line.components.empty()) {
            LC::scalar dr = std::sqrt(cell[0] * cell[0] + cell[1] * cell[1] + cell[2] * cell[2]);
            LC::scalar min_point_dist = 2.0 * _widget.knot_interaction_handle.point_density * dr;

            std::vector<Eigen::Vector3d> distinct_points;
            VKH::buildProcessedKnotGraph(_vortex_line.vertices, min_point_dist, cellf, distinct_points);

            if (!distinct_points.empty()) {
                _processedVortexKnot.emplace_back();
                VortexKnot& knot = _processedVortexKnot.back();

                if (!_vortexKnot.empty()) {
                    knot.knotColor = _vortexKnot.front().knotColor;
                }

                std::function<Vector3(void*, std::size_t)> access = [](void* points, std::size_t i) {
                    const auto* pt = static_cast<const Eigen::Vector3d*>(points) + i;
                    return Vector3(pt->x(), pt->y(), pt->z());
                    };

                knot.points.polyRadius = 0.5 * min_point_dist;
                knot.points.Init(&distinct_points[0].x(), access, distinct_points.size());
            }
        }
        _widget.generateKnots = false;
    }

    LC::Math::Isosurface<float*, float> gen;

    if (_vortexShell) {
        auto shell = std::make_unique<float[]>(sub_grid.vol);
        VKH::computeHelicalShellField(sub_grid, field_nn.get(), _vortexShell->queryValue, _vortexShell->isoLevel, shell.get());

        std::array<float, 4> color2 = { 1.f, 1.f, 1.f, _vortexShell->alpha };
        gen.GenerateSurface(shell.get(), _vortexShell->isoLevel, vNew, cell, color2);

        if (gen.isSurfaceValid()) {
            uint nVert = gen.NumSurfaceVertices();
            uint nInd = gen.NumSurfaceIndices();
            std::unique_ptr<LC::Math::IsoVertex[]> verts(gen.ReleaseSurfaceVertices());
            std::unique_ptr<uint[]> indices(gen.ReleaseSurfaceIndices());

            auto smoothing_result = TaubinSmoothing(
                verts.get(), nVert, indices.get(), nInd, field_nn.get(),
                vNew, { 1.f, 1.f, 1.f, 1.f }, false, false, 0
            );
            auto& vertices = std::get<0>(smoothing_result);
            auto& triangles = std::get<1>(smoothing_result);

            VKH::colorVortexShellVertices(
                vertices, sub_grid, base_grid, chi_field.get(),
                _quaternion_field.get(), _widget.chiColorScheme, _vortexShell->invertNormals
            );

            LC_INFO("Successfully generated surface (verts = {0}, indices = {1})", nVert, nInd);
            registerSurfaceDrawable(
                *_vortexShell, reinterpret_cast<LC::Surface::Vertex*>(vertices.data()),
                vertices.size(), &triangles[0][0], triangles.size() * 3, _vortexShell->noCull
            );
        }
    }

    auto buildQuaternionColoredSurface = [&](auto& target_iso, float* iso_field) {
        VKH::applyDomainInversion(iso_field, base_grid.vol, target_iso->isoLevel, false);

        std::array<float, 4> color2 = { 1.f, 1.f, 1.f, target_iso->alpha };
        gen.GenerateSurface(iso_field, target_iso->isoLevel, data->voxels, cell, color2);

        if (gen.isSurfaceValid()) {
            uint nVert = gen.NumSurfaceVertices();
            uint nInd = gen.NumSurfaceIndices();
            std::unique_ptr<LC::Math::IsoVertex[]> verts(gen.ReleaseSurfaceVertices());
            std::unique_ptr<uint[]> indices(gen.ReleaseSurfaceIndices());

            VKH::colorVerticesByMaxPion(verts.get(), nVert, base_grid, _quaternion_field.get());

            LC_INFO("Successfully generated surface (verts = {0}, indices = {1})", nVert, nInd);
            registerSurfaceDrawable(
                *target_iso, reinterpret_cast<LC::Surface::Vertex*>(verts.get()),
                nVert, indices.get(), nInd, target_iso->noCull
            );
        }
        };

    if (_baryonIsosurface && _baryon_density) {
        auto iso = std::make_unique<float[]>(base_grid.vol);
        for (uint i = 0; i < base_grid.vol; ++i) {
            iso[i] = std::abs(_baryon_density[i] - _baryonIsosurface->query);
        }
        buildQuaternionColoredSurface(_baryonIsosurface, iso.get());
    }

    if (_quaternion_field && _pionPreimage) {
        auto iso = std::make_unique<float[]>(base_grid.vol);
        double px = std::abs(_pionPreimage->pi1);
        double py = std::abs(_pionPreimage->pi2);
        double pz = std::abs(_pionPreimage->pi3);

        for (uint i = 0; i < base_grid.vol; ++i) {
            const Eigen::Quaterniond& q = _quaternion_field[i];
            double dx = std::abs(q.x()) - px;
            double dy = std::abs(q.y()) - py;
            double dz = std::abs(q.z()) - pz;
            iso[i] = static_cast<float>(q.w() * q.w() + dx * dx + dy * dy + dz * dz);
        }
        buildQuaternionColoredSurface(_pionPreimage, iso.get());
    }
}

/*
    Window to allow the generation of cross sectional planes
*/
void Sandbox::handleMultiplaneWindow() {
    if (_widget.multiplane_window) {
        Dataset* data = (Dataset*)(_solver->GetDataPtr());
        ImGui::Begin("Multiplane Tool", &_widget.multiplane_window);
        // Display multiplanes
        _multiplane_manager.cell = data->cell_dims;
        _multiplane_manager.Gui(data, _manipulator.get(), _transparentShader, _transparentDrawables);
        ImGui::End();
    }
}

void Sandbox::handleLehmannWindow() {
    if (_widget.lehmann_tool_window) {
        Dataset* data = (Dataset*)(_solver->GetDataPtr());
        ImGui::Begin("Lehmann Cluster Manipulator", &_widget.lehmann_tool_window);

        // Resets/initializes lehman cluster
        if (ImGui::Button("Initialize Lehmann cluster")) {
            if (_widget.lehmann_reset_helical_bg) {
                _widget.updateImage = true;
                auto helical = LC::Math::Planar(2, 1);
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

                            auto n = helical(pos[0], pos[1], pos[2] - _widget.helicalLayerOffset);

                            for (int d = 0; d < 3; d++) {
                                nn(i, j, k, d) = n[d];
                            }
                        }
            }
            _widget.updateImage = true;
            _lehman_xsection.fill(data->voxels, data->cell_dims, _widget.lehman_upsample);
            uint vol = data->voxels[0] * data->voxels[1] * data->voxels[2];
            for (uint i = 0; i < _lehman_xsection.size(); i++) {
                uint id = _lehman_xsection.index(i);
                data->directors[id] = _lehman_xsection.nodes[i].director.x();
                data->directors[id + vol] = _lehman_xsection.nodes[i].director.y();
                data->directors[id + 2 * vol] = _lehman_xsection.nodes[i].director.z();
            }
        }
        ImGui::SameLine();
        ImGui::Checkbox("Reset Helical field##LehmannWindow", &_widget.lehmann_reset_helical_bg);

        ImGui::PushItemWidth(75.f);
        ImGui::SameLine();
        ImGui::InputInt("Up-sample##Lehman", &_widget.lehman_upsample);

        ImGui::InputFloat("Arclength (dz)", &_widget.lehmanArclength);
        ImGui::SameLine();
        ImGui::InputFloat("z translation (dz)", &_widget.lehman_dz);

        ImGui::InputFloat("Total z distance (p)", &_widget.total_lehman_z_dist);
        ImGui::SameLine();
        ImGui::InputFloat("Total fwd distance (p)", &_widget.total_lehman_forward_dist);
        ImGui::PopItemWidth();

        if (ImGui::Button("Forward")) {
            _widget.updateImage = true;
            uint vol = data->voxels[0] * data->voxels[1] * data->voxels[2];
            LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);
            int iterations = _widget.total_lehman_forward_dist / (0.5 * dz);
            for (int it = 0; it < iterations; it++) {
                // Apply the translation
                _lehman_xsection.forward(0.5 * dz);
                // Fill data
                for (uint i = 0; i < _lehman_xsection.size(); i++) {
                    uint id = _lehman_xsection.index(i);
                    data->directors[id] = _lehman_xsection.nodes[i].director.x();
                    data->directors[id + vol] = _lehman_xsection.nodes[i].director.y();
                    data->directors[id + 2 * vol] = _lehman_xsection.nodes[i].director.z();
                }
            }
        }

        ImGui::SameLine();

        bool z_translation = false;
        if (ImGui::Button("Up")) {
            _widget.updateImage = true;
            z_translation = true;
        }

        ImGui::SameLine();

        int sgn = 1;

        if (ImGui::Button("Down")) {
            _widget.updateImage = true;
            sgn = -1;
            z_translation = true;
        }


        // Fill data with translated lehman cluster
        if (z_translation) {
            uint vol = data->voxels[0] * data->voxels[1] * data->voxels[2];

            LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);
            int iterations = _widget.total_lehman_z_dist / (_widget.lehman_dz * dz);

            for (int it = 0; it < iterations; it++) {
                // Apply the rotation
                _lehman_xsection.rotation(_widget.lehmanArclength * dz, sgn * _widget.lehman_dz * dz);
                // Fill data
                for (uint i = 0; i < _lehman_xsection.size(); i++) {
                    uint id = _lehman_xsection.index(i);
                    data->directors[id] = _lehman_xsection.nodes[i].director.x();
                    data->directors[id + vol] = _lehman_xsection.nodes[i].director.y();
                    data->directors[id + 2 * vol] = _lehman_xsection.nodes[i].director.z();
                }
            }
        }
        ImGui::End();
    }

}

/*
    Export the vortex knots produced as a set of vertices and links (only works on well-defined knots--no multiconnected vertices)
    to be read in KnotPlot.
*/
std::tuple<std::vector<Eigen::Vector3d>, std::vector<int>> exportVortexKnot(
    const std::vector<MeshLib::PNCVertex<float>>& vertices,
    const std::vector<MeshLib::Triangle>& triangles,
    LC::scalar min_point_dist,
    const std::string& fname,
    bool append)
{
    std::vector<Eigen::Vector3d> points;
    std::vector<int> npoints_per_component;
    std::map<uint, Face> unvisited_map;

    const uint nTriangles = static_cast<uint>(triangles.size());
    for (uint tri = 0; tri < nTriangles; ++tri) {
        unvisited_map.insert({ tri, Face(triangles[tri][0], triangles[tri][1], triangles[tri][2]) });
    }

    auto vortex_components = find_all_components_graph(unvisited_map, vertices.size(), 50);
    LC_INFO("Vortex components detected - {0}", vortex_components.size());

    if (vortex_components.empty()) {
        return {};
    }

    std::vector<std::vector<Eigen::Vector3d>> distinct_component_points;
    const LC::scalar minDistSq = min_point_dist * min_point_dist;

    for (std::size_t c = 0; c < vortex_components.size(); ++c) {
        std::vector<Eigen::Vector3d> pointset;
        for (auto id : vortex_components[c]) {
            const auto& v = vertices[id].position;
            Eigen::Vector3d candidate(v[0], v[1], v[2]);

            bool valid_pt = true;
            for (const auto& pt : pointset) {
                if ((candidate - pt).squaredNorm() < minDistSq) {
                    valid_pt = false;
                    break;
                }
            }

            if (valid_pt) {
                pointset.emplace_back(candidate);
            }
        }

        LC_INFO("Component {1}: {0} distinct points detected", pointset.size(), c);
        if (pointset.size() < 5) {
            LC_WARN("Not enough points detected, discarding component");
            continue;
        }
        distinct_component_points.push_back(std::move(pointset));
    }

    constexpr int knn = 5;
    LC::Math::Metric<LC::scalar> metric;
    metric.Bcs = { false, false, false };

    const LC::scalar maxEdgeDist = 2.0 * min_point_dist;
    const LC::scalar cosStraightThreshold = -0.35; // ~110 deg to 180 deg continuation

    for (auto& distinct_points : distinct_component_points) {
        const std::size_t nDistinct = distinct_points.size();
        auto positions = std::make_unique<LC::scalar[]>(nDistinct * 3);
        auto query_domain = std::make_unique<std::size_t[]>(nDistinct);
        auto neighbors = std::make_unique<std::size_t[]>(nDistinct * knn);

        for (std::size_t i = 0; i < nDistinct; ++i) {
            query_domain[i] = i;
            for (int d = 0; d < 3; ++d) {
                positions[i + d * nDistinct] = distinct_points[i](d);
            }
        }

        LC::Algorithm::knn_c(
            positions.get(), nDistinct,
            query_domain.get(), nDistinct,
            metric, knn, static_cast<LC::scalar*>(nullptr), neighbors.get()
        );

        std::vector<Eigen::Vector3d> parsed_points;
        parsed_points.reserve(nDistinct);

        for (std::size_t i = 0; i < nDistinct; ++i) {
            const std::size_t n1 = neighbors[nDistinct + i];
            const LC::scalar dist = (distinct_points[i] - distinct_points[n1]).norm();
            if (dist <= maxEdgeDist) {
                parsed_points.emplace_back(distinct_points[i]);
            }
        }

        const std::size_t nParsed = parsed_points.size();
        if (nParsed < static_cast<std::size_t>(knn)) {
            continue;
        }

        if (nParsed != nDistinct) {
            LC_INFO("Discarding {0} isolated points", nDistinct - nParsed);
            for (std::size_t i = 0; i < nParsed; ++i) {
                query_domain[i] = i;
                for (int d = 0; d < 3; ++d) {
                    positions[i + d * nParsed] = parsed_points[i](d);
                }
            }

            LC::Algorithm::knn_c(
                positions.get(), nParsed,
                query_domain.get(), nParsed,
                metric, knn, static_cast<LC::scalar*>(nullptr), neighbors.get()
            );
        }

        Graph graph(static_cast<int>(nParsed));

        struct TangentCandidate {
            uint idx;
            Eigen::Vector3d dir;
            LC::scalar dist;
        };

        for (std::size_t i = 0; i < nParsed; ++i) {
            std::vector<TangentCandidate> tangents;
            tangents.reserve(knn - 1);

            for (int k = 1; k < knn; ++k) {
                const std::size_t nb_k = neighbors[k * nParsed + i];
                if (nb_k >= nParsed || nb_k == i) continue;

                Eigen::Vector3d diffVec = parsed_points[nb_k] - parsed_points[i];
                const LC::scalar dist = diffVec.norm();

                if (dist > 0.0 && dist < maxEdgeDist) {
                    tangents.push_back({ static_cast<uint>(nb_k), diffVec / dist, dist });
                }
            }

            if (tangents.empty()) {
                continue;
            }
            else if (tangents.size() == 1) {
                graph.addEdge(static_cast<uint>(i), tangents[0].idx);
            }
            else if (tangents.size() == 2) {
                graph.addEdge(static_cast<uint>(i), tangents[0].idx);
                graph.addEdge(static_cast<uint>(i), tangents[1].idx);
            }
            else {
                // Find the pair of neighbors (a, b) that forms the straightest line through i
                // Outward unit vectors from i on the same curve satisfy dir_a . dir_b ~ -1
                int bestA = 0;
                int bestB = 1;
                LC::scalar minDot = (std::numeric_limits<LC::scalar>::max)();

                for (int a = 0; a < static_cast<int>(tangents.size()); ++a) {
                    for (int b = a + 1; b < static_cast<int>(tangents.size()); ++b) {
                        // Weight slightly by distance so closer anti-parallel points win ties
                        const LC::scalar dotVal = tangents[a].dir.dot(tangents[b].dir);
                        if (dotVal < minDot) {
                            minDot = dotVal;
                            bestA = a;
                            bestB = b;
                        }
                    }
                }

                graph.addEdge(static_cast<uint>(i), tangents[bestA].idx);
                graph.addEdge(static_cast<uint>(i), tangents[bestB].idx);

                if (tangents.size() == 4) {
                    std::array<int, 2> rem;
                    int rIdx = 0;
                    for (int c = 0; c < 4; ++c) {
                        if (c != bestA && c != bestB) {
                            rem[rIdx++] = c;
                        }
                    }

                    const LC::scalar secondDot = tangents[rem[0]].dir.dot(tangents[rem[1]].dir);
                    if (secondDot < cosStraightThreshold) {
                        // Bridge the second intersecting strand directly across so the two knots don't merge
                        graph.addEdge(tangents[rem[0]].idx, tangents[rem[1]].idx);
                    }
                }
            }
        }

        auto components = graph.connectedComponents(3);
        for (const auto& comp : components) {
            if (comp.empty()) continue;
            npoints_per_component.push_back(static_cast<int>(comp.size()));
            for (const auto& ci : comp) {
                points.emplace_back(parsed_points[ci]);
            }
        }
    }

    if (!points.empty()) {
        std::ofstream ofile(fname, std::ios::out | std::ios::binary);
        const int numComponents = static_cast<int>(npoints_per_component.size());

        if (ofile.is_open()) {
            ofile.write(reinterpret_cast<const char*>(&numComponents), sizeof(int));
            int offset = 0;
            for (int i = 0; i < numComponents; ++i) {
                ofile.write(reinterpret_cast<const char*>(&npoints_per_component[i]), sizeof(int));
                ofile.write(reinterpret_cast<const char*>(&points[offset]), sizeof(Eigen::Vector3d) * npoints_per_component[i]);
                offset += npoints_per_component[i];
            }
        }

        std::filesystem::path txtPath(fname);
        txtPath.replace_extension(".txt");

        std::ios_base::openmode mode = std::ios::out;
        if (append) mode |= std::ios::app;

        std::ofstream ofile_knotplot(txtPath, mode);
        if (ofile_knotplot.is_open()) {
            int offset = 0;
            for (const int ni : npoints_per_component) {
                for (int i = 0; i < ni; ++i) {
                    const auto& p = points[offset + i];
                    ofile_knotplot << static_cast<float>(p.x()) * 20.0f << " "
                        << static_cast<float>(p.y()) * 20.0f << " "
                        << static_cast<float>(p.z()) * 20.0f << "\n";
                }
                offset += ni;
                ofile_knotplot << "\n";
            }
        }
    }

    return { points, npoints_per_component };
}


/*
    Call total energy function on the solver data
*/
void Sandbox::computeEnergy() {
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    FOFDSolver* solver = (FOFDSolver*)(_solver.get());

    // Update list
    LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
    LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
    LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);
    LC::scalar energy = 0.0;

    if (_widget.radioEn == 0)
        energy = solver->TotalEnergy();
    else if (_widget.radioEn == 1)
        energy = solver->TotalEnergyFunctionalDerivativeAbsSum();

    _widget.energy_series.push_back(energy);
    _widget.series_x_axis.push_back(data->numIterations);

    // Update time series for free energy

    while (_widget.energy_series.size() > _widget.max_graph_points) {
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