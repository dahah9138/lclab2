#include <lclab2.h>
#include "Widget.h"
#include "find_components.h"

using namespace Magnum;
using namespace Math::Literals;
using FOFDSolver = LC::FrankOseen::Electric::FOFDSolver;
using Dataset = FOFDSolver::Dataset;

#define MAX_GRAPH_POINTS 1000
// Current use case for interpolating preimages is inefficient,
// use to interpolate director field THEN copy before passing to isoGenerator
#define TESTINTERP 0
#define USE_RBFINTERP 1

// Smooth isosurfaces
void SmoothIsosurface(LC::Math::IsoVertex* verts, unsigned int* indices, unsigned int nVert, unsigned int nInd, int iterations, float smoothingValue, int smoothingType) {
    unsigned int nTriangles = nInd / 3;

    if (iterations < 1)
        return;

    // Extract vertex data
    Smoothing::Mesh mesh;
    mesh._vertices.resize(nVert);
    mesh._normals.resize(nVert);
    mesh._triangles.resize(nTriangles);

    for (int i = 0; i < nVert; i++) {
        mesh._vertices[i].x = verts[i].position[0];
        mesh._vertices[i].y = verts[i].position[1];
        mesh._vertices[i].z = verts[i].position[2];
    }

    for (int i = 0; i < nTriangles; i++) {
        mesh._triangles[i][0] = indices[3 * i];
        mesh._triangles[i][1] = indices[3 * i + 1];
        mesh._triangles[i][2] = indices[3 * i + 2];
    }

    Smoothing::Vertex_to_1st_ring_vertices first_ring;
    Smoothing::Vertex_to_face v_to_face;
    v_to_face.compute(mesh);
    first_ring.compute(mesh, v_to_face);
    
    if (smoothingType == 0)
        mesh._vertices = Smoothing::explicit_laplacian_smoothing(mesh._vertices, first_ring._rings_per_vertex, iterations, smoothingValue);
    else if (smoothingType == 1)
        mesh._vertices = Smoothing::smooth_iterative(mesh._vertices, first_ring._rings_per_vertex, iterations, smoothingValue);
    else if (smoothingType == 2)
        mesh._vertices = Smoothing::implicit_laplacian_smoothing(mesh._vertices, first_ring._rings_per_vertex, iterations, smoothingValue);

    // Copy data back
    for (int i = 0; i < nVert; i++) {
        verts[i].position[0] = mesh._vertices[i].x;
        verts[i].position[1] = mesh._vertices[i].y;
        verts[i].position[2] = mesh._vertices[i].z;
    }

    // Set normals to zero
    for (unsigned int i = 0; i < nVert; i++) {
        verts[i].normal[0] = 0;
        verts[i].normal[1] = 0;
        verts[i].normal[2] = 0;
    }

    // Recompute normals
    for (int i = 0; i < nTriangles; i++) {

        LC::Math::VECTOR3D vec1, vec2, normal;
        unsigned int id0, id1, id2;
        id0 = indices[i * 3];
        id1 = indices[i * 3 + 1];
        id2 = indices[i * 3 + 2];
        vec1[0] = verts[id1].position[0] - verts[id0].position[0];
        vec1[1] = verts[id1].position[1] - verts[id0].position[1];
        vec1[2] = verts[id1].position[2] - verts[id0].position[2];
        vec2[0] = verts[id2].position[0] - verts[id0].position[0];
        vec2[1] = verts[id2].position[1] - verts[id0].position[1];
        vec2[2] = verts[id2].position[2] - verts[id0].position[2];

        LC::Math::CrossProduct(vec1, vec2, normal);

        verts[id0].normal[0] += normal[0];
        verts[id0].normal[1] += normal[1];
        verts[id0].normal[2] += normal[2];
        verts[id1].normal[0] += normal[0];
        verts[id1].normal[1] += normal[1];
        verts[id1].normal[2] += normal[2];
        verts[id2].normal[0] += normal[0];
        verts[id2].normal[1] += normal[1];
        verts[id2].normal[2] += normal[2];
    }

    // Normalize normals
    for (unsigned int i = 0; i < nVert; i++) {
        float len = 0.0f;
        for (int d = 0; d < 3; d++)
            len += verts[i].normal[d] * verts[i].normal[d];

        len = sqrt(len);

        verts[i].normal[0] /= len;
        verts[i].normal[1] /= len;
        verts[i].normal[2] /= len;
    }
}

enum class Axis { x = 0, y = 1, z = 2 };
struct CrossX {
    Float axisPosition;
    Axis axis;
    std::pair<Containers::Optional<GL::Mesh>, std::unique_ptr<LC::DynamicColorSheet>> section;
    Containers::Optional<LC::EllipsoidArray> nematic;
    bool draw = true;
    bool draw_nematic = true;
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

struct PionPreimage {
    PionPreimage() {}
    PionPreimage(float p1, float p2, float p3, float iso = 0.07f) : pi1(p1), pi2(p2), pi3(p3), isoLevel(iso) {}
    PionPreimage(const std::array<float, 3> &pi, float iso = 0.07f) : pi1(pi[0]), pi2(pi[1]), pi3(pi[2]), isoLevel(iso) {}

    friend bool operator == (const PionPreimage& p1, const PionPreimage& p2) {
        return (p1.pi1 == p2.pi1 && p1.pi2 == p2.pi2 && p1.pi3 == p2.pi3) ? 1 : 0;
    }

    float pi1 = 0.0f;
    float pi2 = 0.0f;
    float pi3 = 1.0f;
    bool draw = true;
    bool noCull = false;
    float isoLevel = 0.07f;
    float alpha = 1.0f;
    bool opentab = true;
    LC::Surface surface;

    Containers::Optional<GL::Mesh> mesh;
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
    std::unique_ptr<short[]> valid_field;
};

struct VortexKnot {
    VortexKnot(float iso = 0.07f) : isoLevel(iso) {}

    friend bool operator == (const VortexKnot& v1, const VortexKnot& v2) {
        return (v1.isoLevel == v2.isoLevel) ? 1 : 0;
    }

    void UpdateColor() {
        for (auto& v : surface.vertices) {
            v.color = knotColor;
        }
        
        mesh = surface.Mesh();
    }


    bool draw = true;
    bool drawSpheres = true;
    float alpha = 1.0f;
    float isoLevel;
    LC::Surface surface;
    Containers::Optional<GL::Mesh> mesh;
    //Containers::Optional<LC::SphereArray> geometry;
    Color4 knotColor{ 1.f, 0.f, 0.f, 1.f };
};

struct VortexShell {
    VortexShell(float iso = 0.07f) : isoLevel(iso) {}

    friend bool operator == (const VortexShell& v1, const VortexShell& v2) {
        return (v1.isoLevel == v2.isoLevel) ? 1 : 0;
    }

    bool draw = true;
    bool noCull = false;
    bool invertNormals = false;
    float alpha = 0.5f;
    float isoLevel;
    float queryValue = 0.0f;
    LC::Surface surface;
    Containers::Optional<GL::Mesh> mesh;
};

struct BaryonIsosurface {

    friend bool operator == (const BaryonIsosurface& v1, const BaryonIsosurface& v2) {
        return (v1.isoLevel == v2.isoLevel) ? 1 : 0;
    }

    bool draw = true;
    bool noCull = false;
    float alpha = 0.5f;
    float isoLevel = 1.999e-2f;
    float query = 2e-2f;
    LC::Surface surface;
    Containers::Optional<GL::Mesh> mesh;
};

struct Zprofile {
    struct Graph {
        std::unique_ptr<Magnum::Vector3[]> data;
        LC::SphereArray grid;
        // Box visual to generate and draw only if z profile window is open
        int vox_x, vox_y;
    };

    struct Box {
        Vector3 dims, translation;
    };

    void GenerateProfile(const LC::scalar* nn, std::array<int,3> vox, std::array<LC::scalar,3> cell, float alpha = 0.0f, int chain_units = 30, float temp = 298.0f, float omegabar = 1.0f, float inversion_temp = 300.0f) {
        
        std::mt19937 generator;
        std::uniform_real_distribution<float> uniform_distribution(-1.0, 1.0);
        auto my_rand = std::bind(uniform_distribution, generator);
        
        graph = Graph{};

        graph->vox_x = vox[0];
        graph->vox_y = vox[1];

        unsigned int plane = vox[0] * vox[1];
        unsigned int vol = plane * vox[2];

        auto flat_idx = [&](int x, int y, int z) {
            x = x >= vox[0] ? x - vox[0] : (x < 0 ? x + vox[0] : x);
            y = y >= vox[1] ? y - vox[1] : (y < 0 ? y + vox[1] : y);
            z = z >= vox[2] ? z - vox[2] : (z < 0 ? z + vox[2] : z);

            return (unsigned int)(x + y * vox[0] + z * plane);
        };

        graph->data = std::unique_ptr<Magnum::Vector3[]>(new Magnum::Vector3[plane]);

        // Amount of thermal randomization
        float tanalpha = tan(alpha);

        float dx = cell[0] / (vox[0] - 1);
        float dy = cell[1] / (vox[1] - 1);
        float dz = cell[2] / (vox[2] - 1);

        // Determine how large the spheres should be:
        graph->grid.polyRadius = 1.f/3.f * sqrt(dx * dx + dy * dy);

        //float dz_par2 = pow(dz * pow(1 - abs(1 - temp / inversion_temp), contraction_factor), 2);
        float dz_perp2 = pow(dz * pow(temp / inversion_temp, expansion_factor), 2);

        float dz_par2 = pow(dz * contraction_factor, 2);
        //float dz_perp2 = pow(dz * expansion_factor, 2);

        // Evaulate each point in the xy plane
        for (int x = 0; x < vox[0]; x++) {
            for (int y = 0; y < vox[1]; y++) {

                // Initialize data to z = -cell[2]/2
                graph->data[x + y * vox[0]] = Magnum::Vector3( -cell[0]/2. + x * dx, -cell[1] / 2. + y * dy, -cell[2]/2.0f );

                float wt = 1.0f;

                for (int z = 0; z < vox[2]; z++) {

                    unsigned int idx = flat_idx(x, y, z);
                    
                    // Find a new director with some thermal noise

                    // Choose a random vector p
                    Eigen::Vector3f p{ 1.f,1.f,1.f };
                    Eigen::Vector3f dir(nn[idx], nn[idx + vol], nn[idx + 2 * vol]);
                    float tolerance = 1e-6;

                    // Make sure it is a valid vector
                    do {
                        // Randomize vector
                        for (int d = 0; d < 3; d++) {
                            
                            p[d] = my_rand();
                        }
                    } while (p.norm() < tolerance || p.cross(dir).norm() < tolerance);

                    Eigen::Vector3f perp = p.cross(dir);
                    perp.normalize();

                    Eigen::Vector3f v = dir + tanalpha * perp;
                    v.normalize();

                    float nz = v[2];

                    float ct2 = nz * nz; // (zhat dot n) = cos(theta) = nz
                    float st2 = 1.0f - ct2;

                    // Not as good as using SOP to define expansion/contraction!
                    float dz_perp2;// = pow(dz * pow(temp / inversion_temp, expansion_factor), 2);
                    float dz_par2; // = pow(dz * contraction_factor, 2);

                    float avgcos2 = 0.0f;
                    for (int m = -1; m < 1; m += 2) {
                        for (int n = -1; n < 1; n += 2) {
                            for (int p = -1; p < 1; p += 2) {

                                unsigned int idx2 = flat_idx(x + m, y + n, z + p);
                                Eigen::Vector3f dir2(nn[idx2], nn[idx2 + vol], nn[idx2 + 2 * vol]);
                                float dotprod = dir.dot(dir2);
                                avgcos2 += dotprod * dotprod / 6.f;
                            }
                        }
                    }

                    // <P2(cos(theta))>
                    float S = 0.5f * (3.0f * avgcos2 - 1.0f);
                    float omega = omegabar * 1e-12 * 4.5e-9; // elastic constant K * (defect core size)
                    double kbT = temp * 1.23e-23;
                    float t = temp / inversion_temp - 1.;
                    float nu = expansion_factor / (1. + exp(-chain_units * t));
   
                    float tuningTheta = 2. * M_PI * dz;
                    float S0 = 0.5f * (3.f * pow(cos(tuningTheta), 2) - 1.f);
           
                    if (abs(S) > S0)
                        S = 1000000;

                    //wt = 1.0 + exp(-omega * pow(S, 2) / kbT);

                    // http://www-f1.ijs.si/~rudi/sola/Seminar-mehka.pdf
                    dz_perp2 = pow(dz * (1. + nu), 2);
                    dz_par2 = pow(dz, 2);
                    

                    float Delta_z = (1.0f - UV_intensity) * dz + UV_intensity * sqrt(ct2 * dz_par2 + st2 * dz_perp2);

                    graph->data[x + y * vox[0]][2] += wt * Delta_z;
                }
                // Reflection/transmission for that column based on interferometry
                //graph->data[x + y * vox[0]][2] *= wt;
            }
        }

        // Now we can initialize the sphere array
        std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* data, std::size_t i) {
            Magnum::Vector3* pos = (Magnum::Vector3*)data;
            return pos[i];
        };

        graph->grid.Init((void*)graph->data.get(), access, plane);
    }
    void Draw(const Magnum::Matrix4 &viewMatrix, const Magnum::Matrix4 &projectionMatrix) {

        graph->grid.Draw(viewMatrix, projectionMatrix);
    }
    Containers::Optional<Graph> graph;
    Box box;
    float expansion_factor = 1.25f; // l_perp
    float contraction_factor = 0.99f; // l_parallel
    float UV_intensity = 1.0f;
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
    void handleZProfileWindow();
    void generatePionTriplet();
    void generateIsosurface();
    void computeEnergy();
    void repeatVolume(int i);

    void findVortexKnotComponents();

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
#if TESTINTERP
    LC::Math::Isosurface<LC::Math::Interp3Map<float>, float> _isoGenerator;
#else
    LC::Math::Isosurface<float*, float> _isoGenerator;
#endif
    std::list<Preimage> _preimages;
    std::list<VortexKnot> _vortexKnot, _processedVortexKnot;
    Containers::Optional<VortexShell> _vortexShell;
    Containers::Optional<PionPreimage> _pionPreimage;
    Containers::Optional<BaryonIsosurface> _baryonIsosurface;
    VortexLine _vortex_line;

    std::unique_ptr<Eigen::Quaternion<LC::scalar>[]> _quaternion_field;
    Containers::Optional<NematicVisual> _quaternion_plane;

    LC::Imaging::ImageSeries _image_series, _nonlin_image_series;
};

Sandbox::Sandbox(const Arguments& arguments) : LC::Application{ arguments,
                                                                Configuration{}.setTitle("FO Electric Simulation")
                                                                               .setWindowFlags(Configuration::WindowFlag::Resizable)

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

            if (ImGui::Button("Quaternion field")) {
                generatePionTriplet();
            }

            if (_quaternion_plane) {
                ImGui::Checkbox("lambda", &_quaternion_plane->drawlambda);
                ImGui::SameLine();
                ImGui::Checkbox("chi", &_quaternion_plane->drawchi);
                ImGui::SameLine();
                ImGui::Checkbox("tau", &_quaternion_plane->drawtau);
            }


            


            // Dropdown menu for lc types
            {
                std::map<LC::FrankOseen::LC_TYPE, std::string> lcMap = LC::FrankOseen::LiquidCrystal::Map();
                dropDownMenu<LC::FrankOseen::LC_TYPE>("LC Type", data->lc_type, lcMap);

                // If custom add customization to gui
                if (data->lc_type == LC::FrankOseen::LC_TYPE::CUSTOM) {
                    ImGui::PushItemWidth(75.f);
                    ImGui::InputDouble("k11", &data->k11.first);
                    ImGui::SameLine();
                    ImGui::InputDouble("k22", &data->k22.first);
                    ImGui::SameLine();
                    ImGui::InputDouble("k33", &data->k33.first);
                    ImGui::InputDouble("epar", &data->epar);
                    ImGui::SameLine();
                    ImGui::InputDouble("eper", &data->eper);
                    ImGui::InputDouble("n0", &data->n0);
                    ImGui::SameLine();
                    ImGui::InputDouble("ne", &data->ne);
                    ImGui::PopItemWidth();
                }
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

            handleZProfileWindow();

            {
                ImGui::SliderFloat("Plane alpha", &_widget.alpha, 0.0f, 1.0f);
                ImGui::Checkbox("Enable POM", &_widget.POM);
                ImGui::SameLine();
                ImGui::Checkbox("Show cell boundary", & _widget.showCellBd);

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

                ImGui::SameLine();

                if (ImGui::Button("Nematic"))
                    ImGui::OpenPopup("nematic_plane_selector");

                ImGui::SameLine();
                ImGui::Checkbox("Coupled", &_widget.couplePlaneAndNematic);


                std::map<std::pair<std::string, Axis>, bool&> planes{ {{"yz", Axis::x }, _crossSections[0].draw },
                    {{"xz", Axis::y }, _crossSections[1].draw },
                    {{"xy", Axis::z }, _crossSections[2].draw } };

                std::map<std::pair<std::string, Axis>, bool&> nematic_planes{ {{"yz", Axis::x }, _crossSections[0].draw_nematic },
                    {{"xz", Axis::y }, _crossSections[1].draw_nematic },
                    {{"xy", Axis::z }, _crossSections[2].draw_nematic } };

                if (ImGui::BeginPopup("plane_selector")) {
                    ImGui::Text("Draw planes");
                    ImGui::Separator();
                    for (auto& p : planes) {
                        ImGui::MenuItem(p.first.first.c_str(), "", &p.second);
                    }
                    ImGui::EndPopup();
                }

                if (ImGui::BeginPopup("nematic_plane_selector")) {
                    ImGui::Text("Draw nematic planes");
                    ImGui::Separator();
                    auto map = planes;
                    if (!_widget.couplePlaneAndNematic) {
                        map = nematic_planes;
                    }
                    for (auto& p : map) {
                        ImGui::MenuItem(p.first.first.c_str(), "", &p.second);
                    }
                    ImGui::EndPopup();
                }

                if (_widget.couplePlaneAndNematic) {
                    for (int i = 0; i < 3; i++)
                        _crossSections[i].draw_nematic = _crossSections[i].draw;
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
                _manipulator->setTransformation(Matrix4::rotationZ(Rad(-M_PI / 2.0f)) * Matrix4::rotationY(Rad(-M_PI / 2.0f)));

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

            //ImGui::Checkbox("GPU", &_widget.GPU);
            //ImGui::SameLine();
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


    updateBoxes();
    drawBoxes();
    


    polyRenderer();

    
    sortObjects(_transparentNormalDrawables);
    if (_widget.drawSurfaces)
        _camera->draw(_transparentNormalDrawables);

    // Sort objects to draw in correct order
    sortObjects(_transparentDrawables);
    _camera->draw(_transparentDrawables);

    // Draw nematic overlay
    for (int id = 0; id < 3; id++) {
        if (_crossSections[id].draw_nematic && _crossSections[id].nematic)
            _crossSections[id].nematic->Draw(_camera->cameraMatrix() * _manipulator->transformationMatrix(), _camera->projectionMatrix());
    }

    // Draw z profile
    if (_widget.drawProfile && _zprofile) {
        _zprofile->Draw(_camera->cameraMatrix()* _manipulator->transformationMatrix(), _camera->projectionMatrix());
    }

    // Draw vortex lines
    //for (auto & line : _vortexKnot)
    //    if (line->geometry && line->drawSpheres)
    //        line->geometry->Draw(_camera->cameraMatrix() * _manipulator->transformationMatrix(), _camera->projectionMatrix());

    if (_quaternion_plane)
        _quaternion_plane->Draw(_camera->cameraMatrix() * _manipulator->transformationMatrix(), _camera->projectionMatrix());

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

template <typename Ty>
Eigen::Quaternion<Ty> fromBasis(Eigen::Matrix<Ty,3,1> a, Eigen::Matrix<Ty, 3, 1> b, Eigen::Matrix<Ty, 3, 1> c) {
    Ty T = a.x() + b.y() + c.z();
    Ty X, Y, Z, W;
    float s;
    if (T > 0) {
        Ty s = sqrt(T + 1) * 2;
        X = (c.y() - b.z()) / s;
        Y = (a.z() - c.x()) / s;
        Z = (b.x() - a.y()) / s;
        W = 0.25 * s;
    }
    else if (a.x() > b.y() && a.x() > c.z()) {
        s = sqrt(1 + a.x() - b.y() - c.z()) * 2;
        X = 0.25 * s;
        Y = (b.x() + a.y()) / s;
        Z = (a.z() + c.x()) / s;
        W = (c.y() - b.z()) / s;
    }
    else if (b.y() > c.z()) {
        s = sqrt(1 + b.y() - a.x() - c.z()) * 2;
        X = (b.x() + a.y()) / s;
        Y = 0.25 * s;
        Z = (c.y() + b.z()) / s;
        W = (b.z() - c.y()) / s;
    }
    else {
        s = sqrt(1 + c.z() - a.x() - b.y()) * 2;
        X = (a.z() + c.x()) / s;
        Y = (c.y() + b.z()) / s;
        Z = 0.25 * s;
        W = (b.x() - a.y()) / s;
    }
    return Eigen::Quaternion<Ty>{W, X, Y, Z};
}

void Sandbox::generatePionTriplet() {
    // Obtain chirality field
    Dataset* data = (Dataset*)_solver->GetDataPtr();

    unsigned int slice = data->voxels[0] * data->voxels[1];
    unsigned int vol = slice * data->voxels[2];

    std::unique_ptr<float[]> nn(new float[3 * vol]);
    for (int i = 0; i < 3 * vol; i++)
        nn[i] = data->directors[i];

    auto vNew = data->voxels;

    unsigned int reduced_slice = vNew[0] * vNew[1];
    unsigned int reduced_vol = reduced_slice * vNew[2];

    // Initialize quaternion field
    using Quatd = Eigen::Quaternion<LC::scalar>;
    _quaternion_field = std::unique_ptr<Quatd[]>(new Quatd[vol]);

    std::unique_ptr<float[]> chi_field;
    std::unique_ptr<short[]> valid_field;

    LC::Math::ChiralityField(nn.get(), chi_field, data->voxels, { (float)data->cell_dims[0],(float)data->cell_dims[1],(float)data->cell_dims[2] }, valid_field, true);
    
    // Extract pion triplet from n and chi
    /*
        Reference triplet: tau    = (1, 0, 0)
                           lambda = (0, 1, 0)
                           chi    = (0, 0, 1)
    */

    for (int i = 0; i < vNew[0]; i++) {
        for (int j = 0; j < vNew[1]; j++) {
            for (int k = 0; k < vNew[2]; k++) {
                
                unsigned int full_idx = i + j * vNew[0] + k * slice;

                Eigen::Vector3f chi, tau;

                Eigen::Vector3f lambda{ nn[full_idx], nn[full_idx + vol], nn[full_idx + 2 * vol] };
                if (!valid_field[full_idx]) {
                    // Make chi perpendicular to n
                    auto zhat = Eigen::Vector3f{ 0.f, 0.f, 1.0f };
                    auto yhat = Eigen::Vector3f{ 0.f, 1.f, 0.f };
                    chi = zhat;
                    
                    auto temp = lambda.cross(chi);
                    if (temp.norm() > 0.)
                        chi = lambda - lambda.dot(zhat) * lambda / lambda.squaredNorm();
                    else {
                        // Guaranteed to be nonzero
                        chi = lambda - lambda.dot(yhat) * lambda / lambda.squaredNorm();
                    }

                    chi.normalize();
                }
                else
                    chi = Eigen::Vector3f{ chi_field[full_idx], chi_field[full_idx + vol], chi_field[full_idx + 2 * vol] };

                tau = -lambda.cross(chi);
                tau.normalize();

                auto q = fromBasis(lambda,tau,chi);

                _quaternion_field[full_idx].x() = q.x();
                _quaternion_field[full_idx].y() = q.y();
                _quaternion_field[full_idx].z() = q.z();
                _quaternion_field[full_idx].w() = q.w();

                // normalize
                _quaternion_field[full_idx].normalize();

                if (q.w() != q.w() || q.x() != q.x() || q.y() != q.y() || q.z() != q.z())
                    LC_WARN("Quaternion NAN");
                
            }
        }
    }

    // Initialize quaternion xy mid-plane
    _quaternion_plane = NematicVisual{};

    // Create plane positions

    float dx = data->cell_dims[0] / (data->voxels[0] - 1);
    float dy = data->cell_dims[1] / (data->voxels[1] - 1);

    // This visualization is not very useful...
    // Need to take the isosurface now
    std::unique_ptr<Magnum::Vector3[]> grid(new Magnum::Vector3[slice]);
    std::unique_ptr<Eigen::Quaternion<float>[]> quaternions(new Eigen::Quaternion<float>[slice]);
    for (int x = 0; x < data->voxels[0]; x++) {
        for (int y = 0; y < data->voxels[1]; y++) {
            grid[x + vNew[0] *y] = Magnum::Vector3(-data->cell_dims[0] / 2. + x * dx, -data->cell_dims[1] / 2. + y * dy, 0.0f);

            auto q = _quaternion_field[x + y * vNew[0] + int(vNew[2] / 2) * slice];

            quaternions[x + vNew[0] * y].x() = q.x();
            quaternions[x + vNew[0] * y].y() = q.y();
            quaternions[x + vNew[0] * y].z() = q.z();
            quaternions[x + vNew[0] * y].w() = q.w();
            
        }
    }

    /*
        Compute the baryon density using the pion triplet
    */

    LC::scalar B = LC::Math::ComputeBaryonDensity(_quaternion_field, _baryon_density, data->voxels);

    LC_INFO("Baryon number = {0}", B);

    /*
        Initialize the midplane of the pion triplet
    */

    _quaternion_plane->dim1 = vNew[0];
    _quaternion_plane->dim2 = vNew[1];

    std::function<Magnum::Vector3(void*, std::size_t)> access = [](void *data, std::size_t i) {
        Magnum::Vector3* mdata = (Magnum::Vector3*)data;
        return mdata[i];
    };

    _quaternion_plane->lambda.Init(grid.get(), access, slice);
    _quaternion_plane->chi.Init(grid.get(), access, slice);
    _quaternion_plane->tau.Init(grid.get(), access, slice);

    // Change colors
    // Use largest magnitude in space components of quaternion field to color according to (r,g,b)
    for (int i = 0; i < slice; i++) {

        
        auto qmat = quaternions[i].toRotationMatrix();

        if (abs(1.0f - quaternions[i].norm()) > 1e-3)
            LC_WARN("Invalid quaternion: norm = {0}", quaternions[i].norm());

        Matrix4 mat;

        for (int ii = 0; ii < 4; ii++) {
            for (int jj = 0; jj < 4; jj++) {
                if (ii == 3 || jj == 3)
                    mat[ii][jj] = 0.f;
                else
                    mat[ii][jj] = qmat(ii, jj);
            }
        }
        mat[3][3] = 1;
                
        float sz = 0.2f * dx;
        Matrix4 transformation;
        {
            transformation = Matrix4::translation(grid[i])
                * mat
                * Matrix4::rotationZ(Rad{ float(M_PI / 2.f) })
                * Matrix4::scaling({ sz,sz,sz });
            _quaternion_plane->lambda.polyInstanceData[i].transformationMatrix = transformation;
            _quaternion_plane->lambda.polyInstanceData[i].normalMatrix = transformation.normalMatrix();
        }

        {
            transformation = Matrix4::translation(grid[i])
                * mat
                * Matrix4::scaling({ sz,sz,sz });
            _quaternion_plane->tau.polyInstanceData[i].transformationMatrix = transformation;
            _quaternion_plane->tau.polyInstanceData[i].normalMatrix = transformation.normalMatrix();
            _quaternion_plane->tau.polyInstanceData[i].color = Color3::red();
        }

        {
            transformation = Matrix4::translation(grid[i])
                * mat
                * Matrix4::rotationX(Rad{ float(M_PI / 2.f) })
                * Matrix4::scaling({ sz,sz,sz });
            _quaternion_plane->chi.polyInstanceData[i].transformationMatrix = transformation;
            _quaternion_plane->chi.polyInstanceData[i].normalMatrix = transformation.normalMatrix();
            _quaternion_plane->chi.polyInstanceData[i].color = Color3::green();
        }
    }


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
        Vector3 scale{ float(data->cell_dims[0]) / (data->voxels[0] - 1),
            float(data->cell_dims[1]) / (data->voxels[1] - 1),
            float(data->cell_dims[2]) / (data->voxels[2] - 1) };

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
                Color3 director_color;

                if (!_widget.nonlinear) director_color = LC::Imaging::Colors::RungeSphere(theta, phi);
                // Green 3 photon nonlinear imaging
                else {
                    if (_widget.nonlinCircular)
                        director_color = Color3::fromHsv({ Deg(120.0f), 1.0f, 1.0f - powf(nz, 6.0f) });
                    else
                        director_color = Color3::fromHsv({ Deg(120.0f), 1.0f, powf(nx * cos(M_PI / 180. * _widget.nonlinTheta) + ny * sin(M_PI / 180. * _widget.nonlinTheta), 6.0f) });
                }
                _crossSections[id].section.second->vertices[cidx].color = { director_color, alpha };
                _crossSections[id].nematic->polyInstanceData[cidx].color = director_color;

                Vector3 translation{ 0.0f, 0.0f, 0.0f };
                translation[id] = data->cell_dims[id] * ((float)hvox / float(data->voxels[id] - 1) - 0.5f);
                translation[xx] = data->cell_dims[xx] / (data->voxels[xx] - 1) * i - data->cell_dims[xx] * .5f;
                translation[yy] = data->cell_dims[yy] / (data->voxels[yy] - 1) * j - data->cell_dims[yy] * .5f;

                _crossSections[id].section.second->vertices[cidx].position[id] = translation[id];
                float hPi = M_PI / 2.f;
                Matrix4 transformation = Matrix4::translation(translation)
                    * Matrix4::rotationZ(Rad{ -hPi + float(phi) })
                    * Matrix4::rotationX(Rad{ hPi - float(theta) })
                    * Matrix4::scaling(0.2f * scale);

                _crossSections[id].nematic->polyInstanceData[cidx].transformationMatrix = transformation;
                _crossSections[id].nematic->polyInstanceData[cidx].normalMatrix = transformation.normalMatrix();

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


    if (data->lc_type != LC::FrankOseen::LC_TYPE::CUSTOM) {
        _pomImager.n0 = LC::FrankOseen::OpticalConstants::LC(data->lc_type, LC::FrankOseen::OpticalConstants::Constant::n_o).first;
        _pomImager.ne = LC::FrankOseen::OpticalConstants::LC(data->lc_type, LC::FrankOseen::OpticalConstants::Constant::n_e).first;
    }
    else {
        _pomImager.n0 = data->n0;
        _pomImager.ne = data->ne;
    }

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
    _widget.radioZProfileIndex = -1;
    
    _zprofile = {};
    _quaternion_plane = {};
    if (_pionPreimage)
        _pionPreimage->draw = false;
    _pionPreimage = {};
    _baryon_density.reset();
    _quaternion_field.reset();

    // Clear components
    _vortex_line.components.clear();

    // Update widget
    for (int d = 0; d < 3; d++) {
        _widget.celldims[d] = data->cell_dims[d];
        _widget.boundaries[d] = data->bc[d];
    }


    // Create a new manipulator and set its parent
    // Reset preimageManipulator as well

    _preimageManipulator = std::make_unique<LC::Drawable::Object3D>();
    _vortexManipulator = std::make_unique<LC::Drawable::Object3D>();
    _processedVortexManipulator = std::make_unique<LC::Drawable::Object3D>();

    _manipulator = std::make_unique<LC::Drawable::Object3D>();
    _manipulator->setParent(&_scene);


    _preimageManipulator->setParent(_manipulator.get());
    _vortexManipulator->setParent(_manipulator.get());
    _processedVortexManipulator->setParent(_manipulator.get());

    // Clear vortexKnot components and processed vortexKnot components
    // Note that since hte manipulators have been reset this does not cause an error
    // to remove components without turning draw off
    _vortexKnot.clear();
    _processedVortexKnot.clear();

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
        if (ax != Axis::z) {
            _crossSections[id].draw = false;
            _crossSections[id].draw_nematic = false;
        }

        _crossSections[id].nematic = LC::EllipsoidArray{};
        unsigned int size = data->voxels[i] * data->voxels[j];
        std::unique_ptr<Vector3[]> positions(new Vector3[size]);
        float dx = data->cell_dims[i] / (data->voxels[i] - 1);
        float dy = data->cell_dims[j] / (data->voxels[j] - 1);

        for (int x = 0; x < data->voxels[i]; x++) {
            for (int y = 0; y < data->voxels[j]; y++) {

                positions[x + y * data->voxels[i]] = { 0.f, 0.f, 0.f };
                positions[x + y * data->voxels[i]][i] = -data->cell_dims[i] / 2. + x * dx;
                positions[x + y * data->voxels[i]][j] = -data->cell_dims[j] / 2. + y * dy;
            }
        }
        // Now we can initialize the sphere array
        std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* data, std::size_t i) {
            Magnum::Vector3* pos = (Magnum::Vector3*)data;
            return pos[i];
        };
        _crossSections[id].nematic->Init(positions.get(), access, size);

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

            if (ImGui::MenuItem("Z-Profile Settings")) {
                _widget.showZProfileWindow = true;
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

        {
            std::array<float, 3> wvlens;
            for (int i = 0; i < 3; i++)
                wvlens[i] = _pomImager.lights[i].wavelength;

            ImGui::InputFloat3("Wavelengths (nm)", &wvlens[0]);

            for (int i = 0; i < 3; i++)
                _pomImager.lights[i].wavelength = wvlens[i];
        }


        {
            std::array<float, 3> intensities;
            for (int i = 0; i < 3; i++)
                intensities[i] = _pomImager.lights[i].intensity;

            ImGui::InputFloat3("Rel. intensity", &intensities[0]);

            for (int i = 0; i < 3; i++)
                _pomImager.lights[i].intensity = intensities[i];
        }

        // Show additional options in a dropdown (TODO)
        ImGui::InputFloat("Gamma", &_pomImager.gamma);
        ImGui::InputFloat("z-depth", &_pomImager.z_scan_ratio);
        ImGui::InputInt("+ layers (p/2)", &_pomImager.additional_layers);

        ImGui::Checkbox("Polarizer", &_pomImager.polarizer);
        ImGui::Checkbox("Analyzer", &_pomImager.analyzer);

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
            ImGui::InputInt("Smoothing iterations", &_widget.smoothingIterations);
            ImGui::SliderFloat("Smoothing value", &_widget.smoothingValue, 0.25f, 200.f);
            ImGui::TextColored({0.f, 1.f, 0.f, 1.f}, "Smoothing type");
            ImGui::RadioButton("Explicit", &_widget.smoothingType, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Iterative", &_widget.smoothingType, 1);
            ImGui::SameLine();
            ImGui::RadioButton("Implicit", &_widget.smoothingType, 2);

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

        ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "------------------------------");

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
            ImGui::RadioButton("Tilt", &_widget.chiColorScheme, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Max quaternion component", &_widget.chiColorScheme, 1);
            ImGui::PopItemWidth();
        }

        if (ImGui::Button("Create Baryon Isosurface")) {
            // Create a new vortex knot
            _baryonIsosurface = BaryonIsosurface{};
        }

        ImGui::SameLine();

        if (ImGui::Button("Remove Baryon Isosurface")) {
            if (_baryonIsosurface)
                _baryonIsosurface->draw = false;
            _baryonIsosurface = {};
        }


        if (_baryonIsosurface) {        
            ImGui::SameLine();
            ImGui::Checkbox("Draw", &_baryonIsosurface->draw);
            ImGui::PushItemWidth(80.0f);
            ImGui::InputFloat("Baryon isovalue", &_baryonIsosurface->isoLevel);
            ImGui::InputFloat("Baryon value", &_baryonIsosurface->query);
            ImGui::PopItemWidth();
        }

        ImGui::PushItemWidth(120.f);

        if (ImGui::Button("Find Components")) {
            findVortexKnotComponents();
        }
          
        ImGui::SameLine();

        if (ImGui::Button("Generate Components")) {
            _widget.generateKnots = true;
        }

        ImGui::PopItemWidth();

        ImGui::PushItemWidth(100.f);
        ImGui::InputInt("Max vortex components", &_widget.maxVortexComponents);
        ImGui::InputInt("Max vortex points", &_vortex_line.maxPoints);
        ImGui::InputInt("Max conic angle (Deg)", &_vortex_line.maxConic);
        ImGui::InputInt("Min component size", &_widget.minComponentSize);
        {
            int knn = _vortex_line.knn;
            ImGui::InputInt("Nearest neighbors", &knn);
            if (knn > 3)
                _vortex_line.knn = knn;
        }
        ImGui::InputFloat("Point separation", &_vortex_line.pointSeparation);
        ImGui::InputFloat("Knot completion threshold", &_widget.knotCompletionDist);
        ImGui::InputInt("Initial cutoff points", &_widget.knotInitialCutoffIterations);
        ImGui::InputInt("Knot refinement iterations", &_widget.knotRefinementIterations);
        ImGui::InputInt("Post-processing up-sampling", &_vortex_line.upSample);
        ImGui::InputFloat("Tube radius", &_vortex_line.tubeRadius);
        ImGui::PopItemWidth();

        {
            int it = 0;

            if (_vortexKnot.size())
                ImGui::TextColored({ 1.f,1.f,0.f,1.f }, "Draw knots");

            for (auto iter = _vortexKnot.begin(); iter != _vortexKnot.end(); iter++) {

                //std::string pts = "Points (K" + std::to_string(it) + ")";
                std::string tub = "Draw (K" + std::to_string(it) + ")";
                std::string col = "K" + std::to_string(it) + " Color";
                std::string rem = "X##K" + std::to_string(it);

                //ImGui::Checkbox(pts.c_str(), &line.drawSpheres);
                //ImGui::SameLine();
                ImGui::Checkbox(tub.c_str(), &iter->draw);
                ImGui::SameLine();
                if (ImGui::ColorEdit4(col.c_str(), &iter->knotColor[0], ImGuiColorEditFlags_NoInputs)) {
                    iter->UpdateColor();
                }

                ImGui::SameLine();

                if (ImGui::Button(rem.c_str())) {
                    iter->draw = false;
                    _vortexKnot.erase(iter);
                    // Remove the corresponding component
                    std::vector<std::vector<uint>>::iterator citer = _vortex_line.components.begin();
                    std::advance(citer, it);
                    _vortex_line.components.erase(citer);
                    break;
                }

                ++it;
            }

            it = 0;

            if (_processedVortexKnot.size())
                ImGui::TextColored({ 1.f,1.f,0.f,1.f }, "Draw processed knots");

            for (auto iter = _processedVortexKnot.begin(); iter != _processedVortexKnot.end(); iter++) {

                //std::string pts = "Points (K" + std::to_string(it) + ")";
                std::string tub = "Draw (K" + std::to_string(it) + ")##processed";
                std::string col = "K" + std::to_string(it) + " Color##processed";
                std::string rem = "X##processedK" + std::to_string(it);
                ++it;

                //ImGui::Checkbox(pts.c_str(), &line.drawSpheres);
                //ImGui::SameLine();
                ImGui::Checkbox(tub.c_str(), &iter->draw);
                ImGui::SameLine();
                if (ImGui::ColorEdit4(col.c_str(), &iter->knotColor[0], ImGuiColorEditFlags_NoInputs)) {
                    iter->UpdateColor();
                }

                ImGui::SameLine();

                if (ImGui::Button(rem.c_str())) {
                    iter->draw = false;
                    _processedVortexKnot.erase(iter);
                    break;
                }

            }
        }

        // Quaternion preimages
        ImGui::PushItemWidth(200.0f);
        ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "Pion field preimages");
        ImGui::SliderFloat3("pi", &_widget.pionComponents[0], 0, 1.0f);
        ImGui::PopItemWidth();

        if (ImGui::Button("Add pion preimage")) {
            _pionPreimage = PionPreimage{_widget.pionComponents, _widget.isoLevel};
        }
        if (_pionPreimage) {
            ImGui::SameLine();
            ImGui::PushItemWidth(50.0f);
            ImGui::SliderFloat("Alpha", &_pionPreimage->alpha, 0.f, 1.f);
            ImGui::PopItemWidth();
            ImGui::SameLine();
            ImGui::Checkbox("Draw", &_pionPreimage->draw);
        }
        

        std::vector <std::array<float, 2>> remove_list;

        for (auto &p : _preimages)
            if (!p.opentab) {
                // Set to not draw (bad things will happen otherwise):
                p.draw = false;
                remove_list.push_back({ p.theta, p.phi });
            }

        // Now we can remove
        for (const auto& p : remove_list) {
            _preimages.remove(Preimage{ p[0], p[1] });
        }
                
        // Check if using global alpha
        if (_widget.global_preimage_alpha)
            for (auto& p : _preimages) 
                p.alpha = _widget.preimage_alpha;


        // Show list to manage preimages
        ImGui::Text("(Theta, Phi) (isolevel)");
        char txtbuffer1[30];
        char txtbuffer2[30];

        if (ImGui::BeginTabBar("Selected preimages", ImGuiTabBarFlags_Reorderable | ImGuiTabBarFlags_AutoSelectNewTabs | ImGuiTabBarFlags_FittingPolicyScroll)) {
            for (auto& p : _preimages) {

                int n1 = sprintf(txtbuffer1, "(%.0f, %.0f) (%f)", p.theta, p.phi, p.isoLevel);
                int n2 = sprintf(txtbuffer2, "(%.0f, %.0f)", p.theta, p.phi);
                // Copy to string
                std::string txt1, txt2;
                for (int i = 0; i < n1; i++) txt1 += txtbuffer1[i];
                for (int i = 0; i < n2; i++) txt2 += txtbuffer2[i];

                if (p.opentab && ImGui::BeginTabItem(txt2.c_str(), &p.opentab, ImGuiTabItemFlags_None))
                {
                    ImGui::Text("%s", txt1.c_str());
                    // Allow isoLevel and alpha to be adjusted
                    ImGui::SliderFloat(("Iso Level " + txt2).c_str(), &p.isoLevel, 0.0f, 0.3f);
                    if (!_widget.global_preimage_alpha)
                        ImGui::SliderFloat(("Iso Alpha " + txt2).c_str(), &p.alpha, 0.0f, 1.0f);

                    ImGui::Checkbox("Draw surface", &p.draw);

                    ImGui::EndTabItem();
                }
            }
            ImGui::EndTabBar();
        }


        if (ImGui::Button("Generate Isosurfaces") || _widget.generateKnots) {
            generateIsosurface();
        }

        ImGui::Checkbox("Draw Surfaces", &_widget.drawSurfaces);

        ImGui::End();

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
                ImPlot::PlotLine(("Z(" + letter + ")").c_str(), ax.get(), line_scan.get(), vox_x);
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
        ImGui::PushItemWidth(75.f);
        ImGui::InputInt3("PBCs", &_widget.boundaries[0]);
        ImGui::InputInt("npp", &_widget.npp);
        ImGui::SameLine();
        ImGui::InputInt("Q", &_widget.topological_charge);
        ImGui::SameLine();
        ImGui::InputInt("Chirality", &_widget.chirality);
        ImGui::PopItemWidth();
        ImGui::RadioButton("Heliknoton", &_widget.hopfion_type, 0);
        ImGui::SameLine();
        ImGui::RadioButton("Hopfion", &_widget.hopfion_type, 1);
        ImGui::SameLine();
        ImGui::RadioButton("Toron", &_widget.hopfion_type, 2);

        int Q = _widget.topological_charge;
        int npp = _widget.npp;
        float v0 = _widget.voltage;

        if (ImGui::Button("Generate")) {

            Dataset::Config cfg;

            if (_widget.hopfion_type == 0) cfg = Dataset::Heliknoton(Q);
            else if (_widget.hopfion_type == 1) cfg = Dataset::Hopfion(Q);
            else cfg = Dataset::Toron();

            data->chirality = _widget.chirality;

            (*data).Voxels(_widget.celldims[0] * npp, _widget.celldims[1] * npp, _widget.celldims[2] * npp)
                .Boundaries(_widget.boundaries[0], _widget.boundaries[1], _widget.boundaries[2])
                .Cell(_widget.celldims[0], _widget.celldims[1], _widget.celldims[2])
                .Configuration(cfg);

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

        if (ImGui::Button("Embed")) {

            // Embed a heliknoton into the selected region
            Dataset::Config cfg;
            if (_widget.hopfion_type == 0) cfg = Dataset::Heliknoton(Q);
            else if (_widget.hopfion_type == 1) cfg = Dataset::Hopfion(Q);
            else cfg = Dataset::Toron();

            (*data).Embed(_widget.shrink_interval_begin, _widget.shrink_interval_end, cfg);

        }

        ImGui::SameLine();


        if (ImGui::Button("Edit existing")) {
            data->chirality = _widget.chirality;
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

        ImGui::PushItemWidth(80.f);
        ImGui::InputInt("Tilt angle (deg)", &_widget.tilt_angle);
        ImGui::SameLine();
        ImGui::InputInt("Tilt direction (deg)", &_widget.tilt_direction);
        ImGui::PopItemWidth();

        if (ImGui::Button("Tilt")) {

            std::unique_ptr<LC::scalar[]> field_nn(new LC::scalar[3 * data->voxels[0] * data->voxels[1] * data->voxels[2]]);

            // tilt direction vector
            LC::scalar tilt_angle_rad = _widget.tilt_angle * M_PI / 180.;
            LC::scalar tilt_dir_rad = _widget.tilt_direction * M_PI / 180.;
            LC::scalar ux = -sin(tilt_dir_rad);
            LC::scalar uy = cos(tilt_dir_rad);
            LC::scalar uz = 0.;
            Eigen::Quaternion<LC::scalar> q(cos(tilt_angle_rad / 2.), ux * sin(tilt_angle_rad / 2.), uy * sin(tilt_angle_rad / 2.), uz * sin(tilt_angle_rad / 2.));
            Eigen::Quaternion<LC::scalar> qbar = q.conjugate();
    
            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
            FOFDSolver::Tensor4 nn_new(field_nn.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

            // Copy
            for (int i = 0; i < data->voxels[0]; i++) {
                for (int j = 0; j < data->voxels[1]; j++) {
                    for (int k = 0; k < data->voxels[2]; k++) {

                        for (int d = 0; d < 3; d++) {
                            nn_new(i, j, k, d) = nn(i, j, k, d);
                        }

                    }
                }
            }


            for (int i = 0; i < data->voxels[0]; i++) {
                for (int j = 0; j < data->voxels[1]; j++) {
                    for (int k = 0; k < data->voxels[2]; k++) {

                        // Rotate using quaternion conjugation
                        Eigen::Quaternion<LC::scalar> n_i(0., nn_new(i, j, k, 0), nn(i, j, k, 1), nn(i, j, k, 2));

                        auto result = q * n_i * qbar;

                        Eigen::Quaternion<LC::scalar> vpos(0., i - (data->voxels[0]-1) / 2., j - (data->voxels[1]-1) / 2., k - (data->voxels[2]-1) / 2.);
                        
                        //auto vresult = q * vpos * qbar;

                        //int in = vresult.x() + (data->voxels[0] - 1) / 2.;
                        //int jn = vresult.y() + (data->voxels[1] - 1) / 2.;
                        //int kn = vresult.z() + (data->voxels[2] - 1) / 2.;

                        int in = i;
                        int jn = j;
                        int kn = k;

                        if (in < 0 || in >= data->voxels[0] || jn < 0 || jn >= data->voxels[1] || kn < 0 || kn >= data->voxels[2]) {
                            continue;
                        }
                        
                        // Normalize result
                        LC::scalar len = sqrt(result.x() * result.x() + result.y() * result.y() + result.z() * result.z());

                        if (len == 0.0) {
                            nn_new(in, jn, kn, 0) = 0.;
                            nn_new(in, jn, kn, 1) = 0.;
                            nn_new(in, jn, kn, 2) = 1.;
                        }
                        else {

                            nn_new(in, jn, kn, 0) = result.x() / len;
                            nn_new(in, jn, kn, 1) = result.y() / len;
                            nn_new(in, jn, kn, 2) = result.z() / len;

                        }

                    }
                }
            }

            // Copy back
            for (int i = 0; i < data->voxels[0]; i++) {
                for (int j = 0; j < data->voxels[1]; j++) {
                    for (int k = 0; k < data->voxels[2]; k++) {
                        
                        for (int d = 0; d < 3; d++) {
                            nn(i, j, k, d) = nn_new(i, j, k, d);
                        }

                    }
                }
            }

            _widget.updateImage = true;
        }

        if (ImGui::Button("Helical field")) {
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

                        auto n = helical(pos[0], pos[1], pos[2]);

                        for (int d = 0; d < 3; d++) {
                            nn(i, j, k, d) = n[d];
                        }
                    }
        }

        ImGui::SameLine();

        if (ImGui::Button("Uniform field")) {
            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

            std::array<float, 3> pos;

            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; i++)
                for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; j++)
                    for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; k++) {
                        nn(i, j, k, 0) = 0;
                        nn(i, j, k, 1) = 0;
                        nn(i, j, k, 2) = 1.0;
                    }
        }

        // Position in [-Lx/2:Lx/2, -Ly/2:Ly/2, -Lz/2:Lz/2]                


        ImGui::TextColored(ImVec4{ 1.0f, 1.0f, 0.0f, 1.0f }, "Translation region");
        ImGui::Text("nn(%d:%d,%d:%d,%d:%d,:)", _widget.shrink_interval_begin[0], _widget.shrink_interval_end[0],
            _widget.shrink_interval_begin[1], _widget.shrink_interval_end[1],
            _widget.shrink_interval_begin[2], _widget.shrink_interval_end[2]
        );
        ImGui::TextColored(ImVec4{ 1.0f, 1.0f, 0.0f, 1.0f }, "Translation direction");
        bool arrow_up = ImGui::ArrowButton("Up", ImGuiDir_Up);
        ImGui::SameLine();
        bool arrow_down = ImGui::ArrowButton("Down", ImGuiDir_Down);
        ImGui::SameLine();
        ImGui::PushItemWidth(100.0f);
        ImGui::InputInt("x", &_widget.translationNumber);
        ImGui::SameLine();
        ImGui::Checkbox("helical", &_widget.helicalTranslation);
        ImGui::PopItemWidth();

        if (arrow_up || arrow_down) {

            FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

            int startz, endz;

            startz = arrow_up ? _widget.shrink_interval_end[2] - 1 : _widget.shrink_interval_begin[2] - 1;
            endz = arrow_up ? _widget.shrink_interval_begin[2] - 1 : _widget.shrink_interval_end[2] - 1;
            int translate = _widget.translationNumber;

            if (arrow_down) translate *= -1;

            auto loop = LC::Utility::create_search_list(LC::Utility::irange_pair(startz, endz));

            LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);
            LC::scalar dtheta = 2. * M_PI * data->chirality * dz;
            LC::scalar cosNtheta = cos(translate * dtheta);
            LC::scalar sinNtheta = sin(translate * dtheta);

            for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; i++) {
                for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; j++) {
                    for (const auto &k : loop) {

                        int kup = k + translate;
                        int kdown = k - translate;

                        kup = kup >= data->voxels[2] ? kup % data->voxels[2] : kup;
                        kup = kup < 0 ? data->voxels[2] + kup : kup;

                        kdown = kdown >= data->voxels[2] ? kdown % data->voxels[2] : kdown;
                        kdown = kdown < 0 ? data->voxels[2] + kdown : kdown;

                        for (int d = 0; d < 3; d++) {

                            // First translate up (down)
                            nn(i, j, kup, d) = nn(i, j, k, d);
                            // then copy from below (above) by same amount
                            nn(i, j, k, d) = nn(i, j, kdown, d);
                        }

                        if (_widget.helicalTranslation) {
                            // Apply a rotation based on chirality, translate, and dz
                            LC::scalar n1x = nn(i, j, kup, 0);
                            LC::scalar n1y = nn(i, j, kup, 1);
                            LC::scalar n2x = nn(i, j, k, 0);
                            LC::scalar n2y = nn(i, j, k, 1);

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

        ImGui::PushItemWidth(100.f);
        ImGui::InputFloat("Separation Distance", &_widget.separationDistance);
        ImGui::InputInt("Theta points", &_widget.interactionThetaPoints);
        ImGui::PopItemWidth();

        if (ImGui::Button("Heliknoton interaction landscape")) {
            // Get all points belonging to the heliknoton
            unsigned int vol = data->voxels[0] * data->voxels[1] * data->voxels[2];

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

            // Normalized handedness to -1
            handedness0(2, 2) = -1.0;

            // Extract the heliknoton
            for (int i = 0; i < data->voxels[0]; i++) {
                for (int j = 0; j < data->voxels[1]; j++) {
                    for (int k = 0; k < data->voxels[2]; k++) {
                        // Diffmag for helical axis
                        unsigned int idx = i + data->voxels[0] * j + data->voxels[0] * data->voxels[1] * k;

                        // Get Chi matrix
                        Eigen::Matrix3d handedness;
                        if ((k <= data->voxels[2] - 3 && k >= 2) && (j <= data->voxels[1] - 3 && j >= 2) && (i <= data->voxels[0] - 3 && i >= 2))
                            handedness = 0.5 / M_PI * LC::Math::HandednessTensor(i, j, k, data->directors.get(), data->voxels, data->cell_dims);
                        else
                            handedness(2, 2) = -1.0;


                        // Compute cost
                        double R2 = 0.;
                        for (int a = 0; a < 3; a++) {
                            for (int b = 0; b < 3; b++) {
                                // R_abR_ab (proportional to "energy" in chi field)
                                R2 += 0.5f * pow(handedness(a, b) - handedness0(a, b), 2);
                            }
                        }
                        
                        // Not part of the background
                        if (R2 > 0.062) {
                            // Add location to the heliknoton
                            Position pos;
                            pos[0] = data->cell_dims[0] * ((LC::scalar)i / (data->voxels[0] - 1) - 0.5);
                            pos[1] = data->cell_dims[1] * ((LC::scalar)j / (data->voxels[1] - 1) - 0.5);
                            pos[2] = data->cell_dims[2] * ((LC::scalar)k / (data->voxels[2] - 1) - 0.5);

                            Director dir;
                            dir[0] = data->directors[idx];
                            dir[1] = data->directors[idx + vol];
                            dir[2] = data->directors[idx + 2 * vol];
                            heliknoton.push_back({ pos, dir });
                        }

                       
                    }
                }
            }

            // Generate the helical background for the cell
            auto helical = LC::Math::Planar(2, 1);
            
            // Create a new volume
            std::array<LC::scalar, 3> cdims{ 8.,8.,8. };
            std::array<int, 3> vNew = { int(_widget.npp * cdims[0]), int(_widget.npp * cdims[1]), int(_widget.npp * cdims[2]) };
            data->directors = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * vNew[0] * vNew[1] * vNew[2]]);
            data->voltage = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);
            data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[vNew[0] * vNew[1] * vNew[2]]);

            FOFDSolver::Tensor4 nn(data->directors.get(), vNew[0], vNew[1], vNew[2], 3);

            // Change voxels and cell dims
            data->voxels = vNew;
            data->cell_dims = cdims;

            auto clear_heliknotons = [&]() {
                // Initialize to the uniform helical bg
                for (int i = 0; i < data->voxels[0]; i++) {
                    for (int j = 0; j < data->voxels[1]; j++) {
                        for (int k = 0; k < data->voxels[2]; k++) {
                            LC::scalar z = data->cell_dims[2] * ((LC::scalar)k / (data->voxels[2] - 1) - 0.5);
                            auto bg = helical(0., 0., z);

                            for (int d = 0; d < 3; d++) {
                                nn(i, j, k, d) = bg[d];
                            }
                        }
                    }
                }
            };

            LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
            LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
            LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);

            auto embed_heliknoton = [&](const Eigen::Vector3d& translation_vector) {
                // Embed the heliknoton in the new volume
                for (const auto &p : heliknoton) {

                    // Apply a rotation + translation to heliknoton position
                    int translatez = translation_vector[2] / dz;
                    LC::scalar dtheta = 2. * M_PI * data->chirality * dz;
                    LC::scalar ctheta = cos(translatez * dtheta);
                    LC::scalar stheta = sin(translatez * dtheta);

                    Position newPos;
                    Director newDir;

                    newPos[0] = ctheta * p.position[0] - stheta * p.position[1];
                    newPos[1] = stheta * p.position[0] + ctheta * p.position[1];
                    newPos[2] = p.position[2];

                    // Rotate director as well
                    newDir[0] = ctheta * p.director[0] - stheta * p.director[1];
                    newDir[1] = stheta * p.director[0] + ctheta * p.director[1];
                    newDir[2] = p.director[2];

                    // Extract new indices
                    int i = (newPos[0] + 0.5 * data->cell_dims[0] + translation_vector[0]) / dx;
                    int j = (newPos[1] + 0.5 * data->cell_dims[1] + translation_vector[1]) / dy;
                    int k = (newPos[2] + 0.5 * data->cell_dims[2] + translation_vector[2]) / dz;

                    // insert a translated+rotated heliknoton
                    for (int d = 0; d < 3; d++)
                        nn(i, j, k, d) = newDir[d];
                }
            };

            // Separation distance
            LC::scalar separation_dist = 1.8;
            std::vector<Eigen::Vector3d> translations;

            struct InteractionPotential {
                LC::scalar theta;
                LC::scalar phi;
                LC::scalar energy;
            };

            struct RadialInteractionPotential {
                LC::scalar radius;
                LC::scalar theta;
                LC::scalar energy;
            };

            // This is the information that will be extracted
            //std::vector<InteractionPotential> interaction_data;
            std::vector<RadialInteractionPotential> interaction_data;
            int thetaPoints = _widget.interactionThetaPoints;
            //int phiPoints = 10;
            int rPoints = 1;

            for (int ti = 0; ti < thetaPoints; ti++) {
                for (int r = 0; r < rPoints; r++) {

                    LC::scalar theta;
                    if (thetaPoints > 1)
                        theta = (LC::scalar)ti / (thetaPoints - 1) * M_PI * 0.5;
                    else
                        theta = M_PI * 0.5;
                    LC::scalar x = (LC::scalar)r / rPoints;
                    LC::scalar rr = _widget.separationDistance;

                    RadialInteractionPotential interaction;
                    interaction.theta = theta;
                    interaction.radius = rr;

                    Eigen::Vector3d displacement;
                    displacement[0] = 0.5 * rr * cos(theta);
                    displacement[1] = 0.5 * rr * sin(theta);
                    displacement[2] = 0.0;

                    translations.emplace_back(displacement);
                    interaction_data.push_back(interaction);
                }
            }

            // Iterate through translations, relaxing and computing energies

            // Get the zero point energy
            LC::scalar density = data->cell_dims[0] * data->cell_dims[1] * data->cell_dims[2] /
                ((data->voxels[0] - 1)* (data->voxels[1] - 1)* (data->voxels[2] - 1));

            clear_heliknotons();
            LC::scalar zero_point_en = solver->TotalEnergy() * density;

            LC_INFO("Zero point energy = {0}", zero_point_en);

#if 0
            int count = 0;
            for (const auto& translation : translations) {
                clear_heliknotons();
                embed_heliknoton(translation);
                embed_heliknoton(-1.0 * translation);
                // Set voltage to 3 V and update voltage 100 times
                solver->SetVoltage(3, 100);
                // Relax for 100 iterations
                solver->Relax(100, true);
                // Compute and save the energy
                interaction_data[count].energy = solver->TotalEnergy() * density - zero_point_en;
                LC_INFO("Energy: {2}, Iteration {0}/{1}", count+1, thetaPoints * rPoints, interaction_data[count].energy);
                ++count;
            }

            // Save data
            std::string saveName = "D:\\dev\\lclab2\\data\\interaction_data_theta_radius_close.bin";
            
            std::ofstream ofile(saveName.c_str(), std::ios::out | std::ios::binary);

            if (!ofile) {
                LC_WARN("Failed to save file <{0}>", saveName.c_str());
            }
            else {
                // Write file version
                ofile.write((char*)&interaction_data[0], interaction_data.size() * sizeof(InteractionPotential));
            }
#else
            // See what the result is without doing the computation
            for (const auto& translation : translations) {
                embed_heliknoton(translation);
                embed_heliknoton(-1.0 * translation);
            }
#endif


            initVisuals();
            _widget.updateImage = true;
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
        repeatVolume(0);
        ImGui::SameLine();
        // Duplicate the volume along the y axis
        repeatVolume(1);
        ImGui::SameLine();
        // Duplicate the volume along the z axis
        repeatVolume(2);



        ImGui::End();
    }

}

void Sandbox::findVortexKnotComponents() {
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    std::array<float, 3> cellf = { data->cell_dims[0], data->cell_dims[1], data->cell_dims[2] };
    std::array<LC::scalar, 3> dr;
    for (int i = 0; i < 3; i++)
        dr[i] = cellf[i] / (data->voxels[i] - 1);

    uint slice = data->voxels[0] * data->voxels[1];
    uint vol = data->voxels[2] * slice;

    // Clean knots
    if (!_vortexKnot.empty()) {
        VortexKnot* knot;
        for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++) {
            knot = &(*it);
            knot->draw = false;
            knot->drawSpheres = false;
        }

        while (!_vortexKnot.empty()) {
            std::list<VortexKnot>::iterator iter;
            for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++)
                iter = it;
            _vortexKnot.erase(iter);
        }
    }

    // Clean manipulator
    _vortexManipulator = std::make_unique<LC::Drawable::Object3D>();
    _vortexManipulator->setParent(_manipulator.get());

    std::unique_ptr<float[]> chi_field, field_nn(new float[3 * vol]);

    for (int i = 0; i < data->voxels[0]; i++) {
        for (int j = 0; j < data->voxels[1]; j++) {
            for (int k = 0; k < data->voxels[2]; k++) {
                uint id = i + data->voxels[0] * j + slice * k;
                field_nn[id] = data->directors[id];
                field_nn[id + vol] = data->directors[id + vol];
                field_nn[id + 2 * vol] = data->directors[id + 2 * vol];
            }
        }
    }

    _vortex_line.numNodes = LC::Math::ChiralityField(field_nn.get(), chi_field, data->voxels, cellf, _vortex_line.valid_field, true, _widget.isoLevel, true);

    std::map<uint, Face> unvisited_list, unvisited_list_local;
    // Key: vertex, value: Triangle
    std::multimap<uint, uint> mmap;

    // Make the initial color white to distinguish components that were too small
    std::array<float, 4> col = { 1.f, 1.f, 1.f, 1.f };

    LC::Math::Isosurface<short*, short> gen2;
    gen2.GenerateSurface(_vortex_line.valid_field.get(), 0, data->voxels, dr, col);

    if (gen2.isSurfaceValid()) {

        unsigned int nInd = gen2.NumSurfaceIndices();
        unsigned int nTriangles = nInd / 3;
        unsigned int nVert = gen2.NumSurfaceVertices();
        LC::Math::IsoVertex* verts = gen2.ReleaseSurfaceVertices();
        unsigned int* indices = gen2.ReleaseSurfaceIndices();

        // Note that verts currently refers to verts and not global indices
        // Therefore, convert indices data to global index space through vertex positions
        for (int tri = 0; tri < nTriangles; tri++) {
            uint i1 = indices[3 * tri];
            uint i2 = indices[3 * tri + 1];
            uint i3 = indices[3 * tri + 2];

            mmap.insert({ i1, tri });
            mmap.insert({ i2, tri });
            mmap.insert({ i3, tri });

            LC::Math::IsoVertex v1 = verts[i1];
            LC::Math::IsoVertex v2 = verts[i2];
            LC::Math::IsoVertex v3 = verts[i3];

            int glob_i1_x = (v1.position[0] / cellf[0] + 0.5) * (data->voxels[0] - 1);
            int glob_i1_y = (v1.position[1] / cellf[1] + 0.5) * (data->voxels[1] - 1);
            int glob_i1_z = (v1.position[2] / cellf[2] + 0.5) * (data->voxels[2] - 1);

            int glob_i2_x = (v2.position[0] / cellf[0] + 0.5) * (data->voxels[0] - 1);
            int glob_i2_y = (v2.position[1] / cellf[1] + 0.5) * (data->voxels[1] - 1);
            int glob_i2_z = (v2.position[2] / cellf[2] + 0.5) * (data->voxels[2] - 1);

            int glob_i3_x = (v3.position[0] / cellf[0] + 0.5) * (data->voxels[0] - 1);
            int glob_i3_y = (v3.position[1] / cellf[1] + 0.5) * (data->voxels[1] - 1);
            int glob_i3_z = (v3.position[2] / cellf[2] + 0.5) * (data->voxels[2] - 1);

            uint glob_i1 = glob_i1_x + glob_i1_y * data->voxels[0] + glob_i1_z * slice;
            uint glob_i2 = glob_i2_x + glob_i2_y * data->voxels[0] + glob_i2_z * slice;
            uint glob_i3 = glob_i3_x + glob_i3_y * data->voxels[0] + glob_i3_z * slice;

            // Store in map
            unvisited_list.insert({ tri, Face(glob_i1, glob_i2, glob_i3) });
            unvisited_list_local.insert({ tri, Face(i1, i2, i3) });
        }

#if 1
        // Find components through the vertices
        _vortex_line.components = find_all_components_graph(unvisited_list_local, nVert, _widget.minComponentSize);

        // Color vertices by their component
        int col_i = 0;
        float Ddeg = 360.f / (float)_vortex_line.components.size();

        // Split apart vertex data
        for (auto& component : _vortex_line.components) {

            Color4 color = Color3::fromHsv({ Deg(Ddeg * col_i++), 1.f, 1.f });
            _vortexKnot.push_back(VortexKnot{});
            VortexKnot* knot;
            {
                std::list<VortexKnot>::iterator iter;
                for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++)
                    iter = it;
                knot = &(*iter);
            }

            knot->knotColor = color;

            // Fill vertex data from component
            std::map<uint, LC::Surface::Vertex> vertices;
            std::vector<uint> indlist;

            for (auto& ci : component) {
                LC::Surface::Vertex v;
                for (int k = 0; k < 3; k++) {
                    v.position[k] = verts[ci].position[k];
                    v.normal[k] = verts[ci].normal[k];
                }
                // Set color
                for (int k = 0; k < 4; k++)
                    v.color[k] = color[k];

                vertices.insert({ ci, v });
            }

            // Parse vertices
            std::vector<LC::Surface::Vertex> vertlist;
            vertlist.reserve(vertices.size());
            for (auto it = vertices.begin(); it != vertices.end(); it++)
                vertlist.emplace_back(it->second);

            std::map<uint, uint> unique_triangles;

            // Now find unique triangles
            for (const auto& vertex : component) {

                // I need a list of triangles corresponding to this vertex
                typedef std::multimap<uint, uint>::iterator MMAPIterator;
                std::pair<MMAPIterator, MMAPIterator> result = mmap.equal_range(vertex);

                // These are all triangles
                for (MMAPIterator it = result.first; it != result.second; it++) {
                    uint tri = it->second;
                    unique_triangles.insert({ tri, tri });
                }
                
            }

            // With unique triangles, now reorder indices
            indlist.reserve(3 * unique_triangles.size());

            for (const auto& face : unique_triangles) {
                uint tri = face.first;
                for (int d = 0; d < 3; d++) {
                    // Vertex index
                    uint ii = indices[3 * tri + d];

                    uint inew = 0;
                    for (auto it = vertices.begin(); it != vertices.end(); it++) {
                        // Found the new location
                        if (it->first == ii) {
                            indlist.push_back(inew);
                            break;
                        }
                        inew++;
                    }
                }
            }

            // Initialize the mesh
            knot->surface.Init(&vertlist[0], vertlist.size(), &indlist[0], indlist.size(), _widget.preimage_translate);

            // Create the mesh
            knot->mesh = knot->surface.Mesh();

            // Add default vortex line to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_vortexManipulator, _phongShader, *(knot->mesh), knot->draw, _transparentNormalDrawables };

            // Convert component indices to global indices

            for (uint i = 0; i < component.size(); i++) {
                LC::Math::IsoVertex v = verts[i];

                int glob_x = (v.position[0] / cellf[0] + 0.5) * (data->voxels[0] - 1);
                int glob_y = (v.position[1] / cellf[1] + 0.5) * (data->voxels[1] - 1);
                int glob_z = (v.position[2] / cellf[2] + 0.5) * (data->voxels[2] - 1);

                component[i] = glob_x + glob_y * data->voxels[0] + glob_z * slice;
            }

        }



#else

        // Find the components
        _vortex_line.components = find_all_components(unvisited_list, _widget.minComponentSize);

        // Color vertices by their component
        int ci = 0;
        float Ddeg = 360.f / (float)_vortex_line.components.size();

        // Split the knot into components
        for (const auto& component : _vortex_line.components) {
            Color4 color = Color3::fromHsv({ Deg(Ddeg * ci++), 1.f, 1.f });
            _vortexKnot.push_back(VortexKnot{});
            VortexKnot* knot;
            {
                std::list<VortexKnot>::iterator iter;
                for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++)
                    iter = it;
                knot = &(*iter);
            }

            knot->knotColor = color;


            // Create component vertices and indices
            std::map<uint, LC::Surface::Vertex> vertices;
            std::vector<uint> indlist;
            indlist.reserve(component.size() * 3);
            for (const auto& face : component) {
                uint tri = face.first;
                for (int d = 0; d < 3; d++) {
                    // Vertex index
                    uint ii = indices[3 * tri + d];
                    LC::Surface::Vertex v;

                    for (int k = 0; k < 3; k++) {
                        v.position[k] = verts[ii].position[k];
                        v.normal[k] = verts[ii].normal[k];
                    }
                    // Set color
                    for (int k = 0; k < 4; k++)
                        v.color[k] = color[k];

                    // Insert if unique
                    vertices.insert({ ii, v });
                }
            }

            // Parse vertices
            std::vector<LC::Surface::Vertex> vertlist;
            vertlist.reserve(vertices.size());
            for (auto it = vertices.begin(); it != vertices.end(); it++)
                vertlist.emplace_back(it->second);

            // With vertlist filled, now find permuted indices
            for (const auto& face : component) {
                uint tri = face.first;
                for (int d = 0; d < 3; d++) {
                    // Vertex index
                    uint ii = indices[3 * tri + d];

                    uint inew = 0;
                    for (auto it = vertices.begin(); it != vertices.end(); it++) {
                        // Found the new location
                        if (it->first == ii) {
                            indlist.push_back(inew);
                            break;
                        }
                        inew++;
                    }
                }
            }

            // Initialize the mesh
            knot->surface.Init(&vertlist[0], vertlist.size(), &indlist[0], indlist.size(), _widget.preimage_translate);

            // Create the mesh
            knot->mesh = knot->surface.Mesh();

            // Add default vortex line to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_vortexManipulator, _phongShader, *(knot->mesh), knot->draw, _transparentNormalDrawables };
        }

#endif

        // Delete data
        delete[] verts;
        delete[] indices;

        LC_INFO("Discovered {0} component(s)", _vortex_line.components.size());
        for (auto& c : _vortex_line.components)
            LC_INFO("- Component size = {0}", c.size());

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



    unsigned int slc = vNew[0] * vNew[1];
    unsigned int size = slc * vNew[2];
    std::unique_ptr<float[]> field_nn(new float[3 * size]);
    FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);
    Eigen::TensorMap<Eigen::Tensor<float, 4>> nn_new(field_nn.get(), vNew[0], vNew[1], vNew[2], 3);
    
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

        // Count the number of points within the preimage
        unsigned int numPointsFound = 0;
        for (unsigned int i = 0; i < size; i++) if (field[i] < pimage.isoLevel) ++numPointsFound;

        // If points found is greater than half of volume, then invert domain
        if (numPointsFound > size / 2) {
            LC_INFO("Preimage points found [{0}/{1}] exceeds half of volume: Inverting domain", numPointsFound, size);
            for (unsigned int i = 0; i < size; i++) {
                if (field[i] > pimage.isoLevel) field[i] = 0.0f;
                else field[i] = 10.0f;
            }
        }


        std::array<float, 4> color;
        {
            auto col = LC::Imaging::Colors::RungeSphere(theta, phi);
            for (int d = 0; d < 3; d++)
                color[d] = col[d];

            color[3] = pimage.alpha;
        }


        _isoGenerator.GenerateSurface(field.get(), pimage.isoLevel, vNew, cell, color);

        if (_isoGenerator.isSurfaceValid()) {

            unsigned int nVert = _isoGenerator.NumSurfaceVertices();
            unsigned int nInd = _isoGenerator.NumSurfaceIndices();

            LC::Math::IsoVertex* verts = _isoGenerator.ReleaseSurfaceVertices();
            unsigned int* indices = _isoGenerator.ReleaseSurfaceIndices();

            SmoothIsosurface(verts, indices, nVert, nInd, _widget.smoothingIterations, _widget.smoothingValue, _widget.smoothingType);

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

    // Make Vortex Knot isosurface
    std::array<float, 4> color1 = { 1., 1., 0., 1. }, color2 = { 1., 1., 1., 1. };
    std::unique_ptr<float[]> chi_field, shell;
    std::unique_ptr<short[]> valid_field;
    LC::Math::Isosurface<float*, float> gen;
    LC::Math::Isosurface<short*, short> gen2;


    unsigned int ill_defined_chi = 0;
    std::array<float,3> cellf = {(float)data->cell_dims[0],(float)data->cell_dims[1],(float)data->cell_dims[2]};

    if (!_widget.chiColorScheme) {
        ill_defined_chi = LC::Math::ChiralityField(field_nn.get(), chi_field, vNew, cellf, valid_field, true, _widget.isoLevel);
    }

    if (_widget.generateKnots) {
        // Clean up all knots excluding components
        if (!_processedVortexKnot.empty()) {
            VortexKnot* knot;
            for (auto it = _processedVortexKnot.begin(); it != _processedVortexKnot.end(); it++) {
                knot = &(*it);
                knot->draw = false;
                knot->drawSpheres = false;
            }

            while (!_processedVortexKnot.empty()) {
                std::list<VortexKnot>::iterator iter;
                for (auto it = _processedVortexKnot.begin(); it != _processedVortexKnot.end(); it++)
                    iter = it;
                _processedVortexKnot.erase(iter);
            }
        }

        // Clean processed knot manipulator
        _processedVortexManipulator = std::make_unique<LC::Drawable::Object3D>();
        _processedVortexManipulator->setParent(_manipulator.get());

        // Define the rbf
        LC::Math::poly_spline<LC::scalar> rbf;
        LC::Math::Metric<LC::scalar> metric;
        metric.SetBCS(LC::Math::PBC({ false, false, false }));
        metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[1]);
        auto vox = data->voxels;
        std::size_t slice = vox[0] * vox[1];
        std::size_t vol = slice * vox[2];


        if (_vortex_line.numNodes && _vortex_line.components.size() > 0) {
            std::vector<Eigen::Vector3d> components_COM;

            // Fill data
            _vortex_line.neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[_vortex_line.knn * _vortex_line.numNodes]);
            _vortex_line.queryDomain = std::unique_ptr<std::size_t[]>(new std::size_t[_vortex_line.numNodes]);
            _vortex_line.positions = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * _vortex_line.numNodes]);

            // key: global index, value: local index
            std::map<unsigned int, unsigned int> queryMap;

            std::unique_ptr<LC::scalar[]> pos_data(new LC::scalar[3 * vol]);
            std::unique_ptr<LC::scalar[]> chi_discriminant(new LC::scalar[vol]);
            
            int ill_idx = 0;
            LC::scalar dx = data->cell_dims[0] / (vox[0] - 1);
            LC::scalar dy = data->cell_dims[1] / (vox[1] - 1);
            LC::scalar dz = data->cell_dims[2] / (vox[2] - 1);

            Eigen::Matrix3d handedness;
            Eigen::Vector3d startpoint;
            unsigned int within_tolerance = 0;

            for (int i = 0; i < vox[0]; i++) {
                for (int j = 0; j < vox[1]; j++) {
                    for (int k = 0; k < vox[2]; k++) {

                        unsigned int sub_idx = i + j * vox[0] + k * slice;

                        // These positions are not correctly translated to the center if there is a cropping...
                        pos_data[sub_idx] = i * dx - data->cell_dims[0] * 0.5;
                        pos_data[sub_idx + vol] = j * dy - data->cell_dims[1] * 0.5;
                        pos_data[sub_idx + 2 * vol] = k * dz - data->cell_dims[2] * 0.5;

                        handedness = 0.5 / M_PI * LC::Math::HandednessTensor(i, j, k, field_nn.get(), vox, cellf);

                        LC::scalar tr2A = pow(handedness.trace(), 2);
                        LC::scalar trA2 = (handedness * handedness).trace();

                        chi_discriminant[sub_idx] = 2. * trA2 - tr2A;

                        // Found a place to interp about with rbf
                        if (!_vortex_line.valid_field[sub_idx]) {
                            _vortex_line.positions[ill_idx] = pos_data[sub_idx];
                            _vortex_line.positions[ill_idx + _vortex_line.numNodes] = pos_data[sub_idx + vol];
                            _vortex_line.positions[ill_idx + 2 * _vortex_line.numNodes] = pos_data[sub_idx + 2 * vol];
                            queryMap.insert({ sub_idx, ill_idx });
                            _vortex_line.queryDomain[ill_idx++] = sub_idx;
                        }
                        else if (abs(1. - chi_discriminant[sub_idx]) < _widget.isoLevel) {
                            ++within_tolerance;
                        }

                        
                    }
                }
            }

            LC_INFO("{0} points within tolerance of {1}% of q^2", within_tolerance, 100.f * _widget.isoLevel);

            for (int loop = 0; loop < min(_vortex_line.components.size(), _widget.maxVortexComponents); loop++) {

                bool loopComponentClosed = false;

                // Choose vertex with the smallest chi discriminant
                uint startId = _vortex_line.components[loop][0];
                LC::scalar smallest_chi = chi_discriminant[startId];

                for (const auto& vertex : _vortex_line.components[loop]) {
                    if (chi_discriminant[vertex] < smallest_chi) {
                        smallest_chi = chi_discriminant[vertex];
                        startId = vertex;
                    }
                }

                startpoint(0) = pos_data[startId];
                startpoint(1) = pos_data[startId + vol];
                startpoint(2) = pos_data[startId + 2 * vol];

                std::unique_ptr<std::size_t[]> nbs(new std::size_t[_vortex_line.numNodes]);

                auto EvaluatePoint = [&](Eigen::Vector3d& p, std::unique_ptr<LC::scalar[]> & pos_subset, std::unique_ptr<LC::scalar[]> &chi_subset, uint sz) {
                    // Find nearest knn neighbors in global grid corresponding to query point
                    LC::Algorithm::knn_c(pos_subset.get(), sz, &p(0), 1, metric, _vortex_line.knn, (LC::scalar*)0, nbs.get());

                    // Position and nearest neighbor data has been filled!
                    // Compute the factorization!
                    _vortex_line.interpolant.ComputeFactorization(pos_subset.get(), nbs.get(), rbf, metric, 1, sz, _vortex_line.knn);

                    // Compute the weights
                    _vortex_line.interpolant.ComputeWeights(chi_subset.get(), nbs.get(), 1, sz, _vortex_line.knn);

                    return _vortex_line.interpolant.Evaluate(0, pos_subset.get(), &p(0), nbs.get(),
                        rbf, metric, 1, sz, _vortex_line.knn);
                };

                // Generate singular line

                std::vector<Eigen::Vector3d> qpoints, qpoints_refined;
                // Starting point
                qpoints.emplace_back(startpoint);

                // Admissible angle of deviation
                LC::scalar deviationAngle = _vortex_line.maxConic;
                LC::scalar deviationAngleRad = deviationAngle * M_PI / 180.;

                // Ball radius to choose next node
                LC::scalar dr = (dx + dy + dz) / 3.;
                LC::scalar radius = _vortex_line.pointSeparation * dr;

                auto Convert3DPointToIdx = [&](const Eigen::Vector3d& pt) {
                    // Convert potential point to index space of pos_data
                    int px = (pt.x() / data->cell_dims[0] + 0.5) * (vox[0] - 1);
                    int py = (pt.y() / data->cell_dims[1] + 0.5) * (vox[1] - 1);
                    int pz = (pt.z() / data->cell_dims[2] + 0.5) * (vox[2] - 1);

                    return px + vox[0] * py + slice * pz;
                };

                auto ConvertIdxTo3DPoint = [&](const unsigned int& idx) {
                    // Convert potential point to index space of pos_data
                    int px = idx / slice;
                    int py = (idx - px * slice) / vox[0];
                    int pz = idx - py * vox[0] - px * slice;

                    Eigen::Vector3d pt((px / LC::scalar(vox[0] - 1) - 0.5) * data->cell_dims[0],
                        (py / LC::scalar(vox[1] - 1) - 0.5) * data->cell_dims[1],
                        (pz / LC::scalar(vox[2] - 1) - 0.5) * data->cell_dims[2]
                    );

                    return pt;
                };

                // Choose second node from most negative neighboring node starting at cube radius * (-1,-1,-1) from current pos

                int resolutionX = 20;
                int resolutionY = 20;
                int resProd = resolutionX * resolutionY;

                std::unique_ptr<Eigen::Vector3d[]> solid_angle_pts(new Eigen::Vector3d[resProd]);


                // f(chi) values at each solid angle point
                std::unique_ptr<LC::scalar[]> solid_angle_f(new LC::scalar[resProd]);

                auto NextPosition = [&](const Eigen::Vector3d& e1, const Eigen::Vector3d& e2, const Eigen::Vector3d& e3, LC::scalar deviation, 
                    std::unique_ptr<LC::scalar[]>& pos_subset, std::unique_ptr<LC::scalar[]>& chi_subset, uint sz) {

                    // Parametrize set of points within solid angle that coincide with pos_data
                    Eigen::Vector3d potential_point;

                    for (int j = 0; j < resolutionY; j++) {
                        for (int i = 0; i < resolutionX; i++) {
                            LC::scalar theta = (j + 1) / LC::scalar(resolutionY + 1) * deviation;
                            LC::scalar phi = i / LC::scalar(resolutionX) * 2. * M_PI;
                            potential_point = qpoints[qpoints.size() - 1]
                                + radius * sin(theta) * cos(phi) * e1
                                + radius * sin(theta) * sin(phi) * e2
                                + radius * cos(theta) * e3;

                            solid_angle_pts[i + resolutionX * j] = potential_point;
                            solid_angle_f[i + resolutionX * j] = EvaluatePoint(potential_point, pos_subset, chi_subset, sz)[0];

                        }
                    }

                    LC::scalar fmin = solid_angle_f[0];
                    unsigned int bestIdx = 0;

                    // Look at all solid angle points for the most negative
                    // Can change this to be Monte Carlo if needed
                    for (int i = 1; i < resProd; i++) {

                        if (solid_angle_f[i] < fmin) {
                            fmin = solid_angle_f[i];
                            bestIdx = i;
                        }
                    }

                    return solid_angle_pts[bestIdx];
                };

                {
                    Eigen::Vector3d xhat{ 1., 0., 0. };
                    Eigen::Vector3d yhat{ 0., 1., 0. };
                    Eigen::Vector3d zhat{ 0., 0., 1. };
                    Eigen::Vector3d bestPos = NextPosition(xhat, yhat, zhat, M_PI, pos_data, chi_discriminant, vol);

                    if (qpoints[0] == bestPos)
                        LC_INFO("Error: Circulation direction is zero");

                    qpoints.emplace_back(bestPos);

                    LC_INFO("Point added (Total = {0}/{1})", qpoints.size(), _vortex_line.maxPoints);
                }


                // With two points, a circulation direction can now be established
                Eigen::Vector3d circulationDirection = qpoints[1] - qpoints[0];
                circulationDirection.normalize();
                // Initial circulation direction (needed to close the loop)
                Eigen::Vector3d circulationDirection0 = circulationDirection;
                Eigen::Vector3d perp1, perp2;

                auto FindOrthonormalBasis = [&]() {
                    perp1 = circulationDirection.cross(Eigen::Vector3d(0., 0., 1.));
                    LC::scalar p1norm2 = perp1.squaredNorm();

                    if (p1norm2) {
                        perp1 = perp1 / sqrt(p1norm2);
                    }
                    else {
                        // 0,1,0 is guaranteed nonzero since circ dir is collinear to 0,0,1
                        perp1 = circulationDirection.cross(Eigen::Vector3d(0., 1., 0.));
                        perp1.normalize();
                    }

                    // Compute perp2
                    perp2 = circulationDirection.cross(perp1);
                };

                // Begin path finding algorithm
                // 1. Look for most optimal direction (most negative) within solid angle about circulation direction
                // - Parametrize solid angle by finding two perpendicular axes to circulationDirection
                // - Look for most optimal point
                // 2. Add point to qpoints
                // 3. Update circulation direction
                // 4. Repeat

                resolutionX = 10;
                resolutionY = 6;
                resProd = resolutionX * resolutionY;

                solid_angle_pts = std::unique_ptr<Eigen::Vector3d[]>(new Eigen::Vector3d[resProd]);


                // f(chi) values at each solid angle point
                solid_angle_f = std::unique_ptr<LC::scalar[]>(new LC::scalar[resProd]);

                uint halfCubeSide = max(2. * _widget.separationDistance, pow(_vortex_line.knn * 3. / 4. / M_PI, 1. / 3.));
                uint cubeSide = 2 * halfCubeSide + 1;
                LC_INFO("Cube side = {0} units", cubeSide);
                uint cubeVol = cubeSide * cubeSide * cubeSide;

                std::unique_ptr<LC::scalar[]> pos_subset(new LC::scalar[3 * cubeVol]);
                std::unique_ptr<LC::scalar[]> chi_subset(new LC::scalar[cubeVol]);

                // Repeat until loop has been fully traversed
                int knotRefinementCount = 0;
                bool cutoff = true;
                while (qpoints.size() < _vortex_line.maxPoints && knotRefinementCount <= _widget.knotRefinementIterations) {
                    // Using circulation direction as the zhat axis, define two orthogonal directions
                    FindOrthonormalBasis();

                    Eigen::Vector3d center = qpoints[qpoints.size() - 1] + radius * circulationDirection;

                    int c0x = (center.x() / cellf[0] + 0.5) * (vox[0] - 1);
                    int c0y = (center.y() / cellf[1] + 0.5) * (vox[1] - 1);
                    int c0z = (center.z() / cellf[2] + 0.5) * (vox[2] - 1);

                    // Fill subset data
                    int xlb = c0x - halfCubeSide;
                    int ylb = c0y - halfCubeSide;
                    int zlb = c0z - halfCubeSide;

                    if (xlb < 0) xlb = 0;
                    if (ylb < 0) ylb = 0;
                    if (zlb < 0) zlb = 0;

                    int xub = (c0x + halfCubeSide) % vox[0];
                    int yub = (c0y + halfCubeSide) % vox[1];
                    int zub = (c0z + halfCubeSide) % vox[2];

                    uint count = 0;
                    auto searchX = LC::Utility::create_search_list(LC::Utility::uirange_pair(xlb, xub));
                    auto searchY = LC::Utility::create_search_list(LC::Utility::uirange_pair(ylb, yub));
                    auto searchZ = LC::Utility::create_search_list(LC::Utility::uirange_pair(zlb, zub));

                    uint realSize = searchX.size() * searchY.size() * searchZ.size();

                    if (realSize != cubeVol) {
                        cubeVol = realSize;
                        pos_subset = std::unique_ptr<LC::scalar[]>(new LC::scalar[3 * cubeVol]);
                        chi_subset = std::unique_ptr<LC::scalar[]>(new LC::scalar[cubeVol]);
                    }

                    for (auto xi : searchX) {
                        for (auto yi : searchY) {
                            for (auto zi : searchZ) {
                                uint i = xi + yi * vox[0] + zi * slice;
                                pos_subset[count] = pos_data[i];
                                pos_subset[count + cubeVol] = pos_data[i + vol];
                                pos_subset[count + 2 * cubeVol] = pos_data[i + 2 * vol];
                                chi_subset[count++] = chi_discriminant[i];
                            }
                        }
                    }


                    Eigen::Vector3d pos = NextPosition(perp1, perp2, circulationDirection, deviationAngleRad,
                        pos_subset, chi_subset, cubeVol);

                    // Update circulation direction and add point
                    circulationDirection = pos - qpoints[qpoints.size() - 1];
                    circulationDirection.normalize();

                    if ((pos - qpoints[0]).norm() <= _widget.knotCompletionDist * radius
                        && circulationDirection0.dot(circulationDirection) > 0.
                        && (qpoints[0] - pos).dot(circulationDirection0) > 0.) {
                        // End of the loop, circulation found closed loop!
                        loopComponentClosed = true;
                        // Remove the initial point and add the final point (better estimate)
                        // Note that this could be used to create a passover of the line
                        // with increased accuracy but it may also not converge...
                        qpoints.erase(qpoints.begin());
                        qpoints.emplace_back(pos);

                        break;
                    }
                    else {
                        // Keep adding points like usual
                        qpoints.emplace_back(pos);

                        // Cut off first 15 points to allow vortex line to converge
                        // to a path first
                        if (cutoff && qpoints.size() > 1 + _widget.knotInitialCutoffIterations) {
                            qpoints.erase(qpoints.begin(), qpoints.begin() + _widget.knotInitialCutoffIterations);
                            // Redefine initial circulation direction
                            circulationDirection0 = qpoints[1] - qpoints[0];
                            circulationDirection0.normalize();
                            cutoff = false;
                            LC_INFO("Cut off {0} points", _widget.knotInitialCutoffIterations);
                        }
                        else {
                            if (qpoints.size() == _vortex_line.maxPoints && _widget.knotRefinementIterations) {
                                qpoints.erase(qpoints.begin());
                                LC_INFO("Knot refinements = {0}/{1}", ++knotRefinementCount, _widget.knotRefinementIterations);
                            }
                            else {
                                LC_INFO("Point added (Total = {0}/{1})", qpoints.size(), _vortex_line.maxPoints);
                            }
                        }

                    }
                }

                // Compute COM
                {
                    Eigen::Vector3d COM(0., 0., 0.);
                    LC::scalar wt = 1./(LC::scalar)qpoints.size();
                    for (const auto& p : qpoints) {
                        COM = COM + p;
                    }
                    COM = wt * COM;

                    bool duplicate = false;

                    for (const auto& com : components_COM) {
                        if ((com - COM).norm() < 2 * radius)
                            duplicate = true;
                    }

                    if (!duplicate)
                        components_COM.push_back(COM);
                    else // Was a duplicate point
                        continue;

                }

                // Vortex knot is valid at this point so generate
                _processedVortexKnot.push_back(VortexKnot{});


                // Refine points via Chaikin's method twice

                for (int it = 0; it < _vortex_line.upSample; it++) {

                    qpoints_refined.resize(2 * (qpoints.size() - !loopComponentClosed));

                    for (int i = 0; i < qpoints.size() - !loopComponentClosed; i++) {

                        int ip1 = (i + 1) % qpoints.size();
                        qpoints_refined[2 * i] = 0.75 * qpoints[i] + 0.25 * qpoints[ip1];
                        qpoints_refined[2 * i + 1] = 0.25 * qpoints[i] + 0.75 * qpoints[ip1];
                    }
                    qpoints = qpoints_refined;
                }

                // Export points
                VortexKnot* knot;
                    
                for (auto it = _processedVortexKnot.begin(); it != _processedVortexKnot.end(); it++) {
                    knot = &(*it);
                }
                    
                //knot->geometry = LC::SphereArray{};

                std::function<Vector3(void*, std::size_t)> access = [](void* points, std::size_t i) {
                    Vector3 result;
                    result[0] = ((Eigen::Vector3d*)points + i)->x();
                    result[1] = ((Eigen::Vector3d*)points + i)->y();
                    result[2] = ((Eigen::Vector3d*)points + i)->z();
                    return result;
                };

                {
                    std::list<VortexKnot>::iterator iter = _vortexKnot.begin();
                    std::advance(iter, loop);
                    knot->knotColor = iter->knotColor;
                }

                //knot->geometry->polyRadius = 0.25f / (float)_vortex_line.upSample * radius;
                //knot->geometry->Init(&qpoints_refined[0].x(), access, qpoints_refined.size());

                // Pass points to tubular surface
                LC::TubularSurface tsurface;
                tsurface.curve = qpoints_refined;
                tsurface.radius = _vortex_line.tubeRadius * 0.5f * dr;
                // Note tubular surface is not compiled, vertices and indices are only generated
                tsurface.Init(loopComponentClosed);

                // Set colors
                for (int i = 0; i < tsurface.data.size(); i++)
                    tsurface.data[i].color = knot->knotColor;

                if (tsurface.data.size() && tsurface.indices.size()) {

                    LC::Math::IsoVertex* verts = (LC::Math::IsoVertex*)&tsurface.data[0];
                    unsigned int* indices = &tsurface.indices[0];
                    // Create the usual surface
                    knot->surface.Init((LC::Surface::Vertex*)verts, tsurface.data.size(),
                        indices, tsurface.indices.size(), _widget.preimage_translate);
                    knot->mesh = knot->surface.Mesh();

                    LC_INFO("Added vortex tube to scene");
                    // Add mesh to the scene
                    new LC::Drawable::TransparentNormalDrawable{ *_processedVortexManipulator, _phongShader, *(knot->mesh), knot->draw, _transparentNormalDrawables };
                }

                
            }

        }
        
        _widget.generateKnots = false;

    }

    if (_vortexShell) {

        color2[3] = _vortexShell->alpha;
        unsigned int vol = vNew[0] * vNew[1] * vNew[2];
        shell = std::unique_ptr<float[]>(new float[vol]);
        Eigen::Matrix3d handedness0, handedness;
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                handedness0(a, b) = 0.0;
            }
        }

        // Normalized handedness to -1
        handedness0(2, 2) = -1.0;

        for (int i = 0; i < vNew[0]; i++) {
            for (int j = 0; j < vNew[1]; j++) {
                for (int k = 0; k < vNew[2]; k++) {
                    // Diffmag for helical axis
                    unsigned int idx = i + vNew[0] * j + vNew[0] * vNew[1] * k;
                    
                    if ((k <= vNew[2] - 3 && k >= 2) && (j <= vNew[1] - 3 && j >= 2) && (i <= vNew[0] - 3 && i >= 2))
                        handedness = 0.5 / M_PI * LC::Math::HandednessTensor(i, j, k, field_nn.get(), vNew, cellf);
                    else
                        handedness(2, 2) = -1.;
                   
                    // Compute cost
                    shell[idx] = 0.0f;
                    double R2 = 0.;
                    for (int a = 0; a < 3; a++) {
                        for (int b = 0; b < 3; b++) {
                            // R_abR_ab (proportional to "energy" in chi field)
                           R2 += 0.5f * pow(handedness(a, b) - handedness0(a, b), 2);
                        }
                    }

                    shell[idx] = abs(R2 - _vortexShell->queryValue);

                    //shell[idx] = pow(chi_field[idx + 2 * size] - 1.0f, 2);
                }
            }
        }

        // Count the number of points within the preimage
        unsigned int helicalPointsFound = 0;

        for (unsigned int i = 0; i < vol; i++) if (shell[i] < _vortexShell->isoLevel) ++helicalPointsFound;

        // If points found is greater than half of volume, then invert domain
        if (helicalPointsFound > vol / 2) {
            LC_INFO("Preimage points found [{0}/{1}] exceeds half of volume: Inverting domain", helicalPointsFound, vol);
            for (unsigned int i = 0; i < vol; i++) {
                if (shell[i] > _vortexShell->isoLevel) shell[i] = 0.0f;
                else shell[i] = 10.0f;
            }
        }


        gen.GenerateSurface(shell.get(), _vortexShell->isoLevel, vNew, cell, color2);

        if (gen.isSurfaceValid()) {

            unsigned int nVert = gen.NumSurfaceVertices();
            unsigned int nInd = gen.NumSurfaceIndices();
            unsigned int nTriangles = nInd / 3;

            LC::Math::IsoVertex* verts = gen.ReleaseSurfaceVertices();
            unsigned int* indices = gen.ReleaseSurfaceIndices();

            SmoothIsosurface(verts, indices, nVert, nInd, 1, 200.f, 2);

            // Go through vertices and color according to vortex shell tilt
            // Note needs to include preimage_translate since embedded volumes can be chosen!

            // Reduced quaternion space SU2/Q8
            Eigen::Matrix3d I3, A1, A2, A3;
            I3 = Eigen::Matrix3d::Identity();
            A1 = Eigen::DiagonalMatrix<LC::scalar, 3, 3>(1., -1., -1.);
            A2 = Eigen::DiagonalMatrix<LC::scalar, 3, 3>(-1., 1., -1.);
            A3 = Eigen::DiagonalMatrix<LC::scalar, 3, 3>(-1., -1., 1.);

            for (int i = 0; i < nVert; i++) {
                // extract indices from position
                int xx = (verts[i].position[0] / data->cell_dims[0] + 0.5) * (vNew[0] - 1);
                int yy = (verts[i].position[1] / data->cell_dims[1] + 0.5) * (vNew[1] - 1);
                int zz = (verts[i].position[2] / data->cell_dims[2] + 0.5) * (vNew[2] - 1);

                // Get chi
                Color3 tilt_color;

                // Chosen color scheme
                // 0 == tilt coloring
                // 1 == quaternion coloring
                if (!_widget.chiColorScheme) {
                    unsigned int idx = xx + yy * vNew[0] + zz * slc;

                    float chi_x = chi_field[idx];
                    float chi_y = chi_field[idx + vol];
                    float chi_z = chi_field[idx + 2 * vol];


                    // Get angle in xy plane
                    {
                        float theta = acos(chi_z);
                        float phi = 2 * atan2(chi_y, chi_x);
                        if (phi < 0.0f) phi = phi + 2 * M_PI;

                        tilt_color = LC::Imaging::Colors::RungeSphere(M_PI / 2., phi);
                        //tilt_color = pow(cos(angle2), 2) * purple + pow(sin(angle2),2) * tilt_color;
                    }
                }
                else {
                    unsigned int idx = xx + yy * data->voxels[0] + zz * data->voxels[0] * data->voxels[1];
                    tilt_color = Color3::blue();
                    Eigen::Quaterniond q{ 0.,0.,0.,1. };
                    if (_quaternion_field)
                        q = _quaternion_field[idx];

                    // Project q onto SU2/Q8
                    auto projection = q.w() * I3 + q.x() * A1 + q.y() * A2 + q.z() * A3;

                    // Make comparisons based on different components

                    // 2q.x() - 2q.y() < 0 ? then y > x
                    int comp = 0;
                    if (projection(0,0) - projection(1,1) > 0.0) {
                        comp = 1;
                    }
                    // 2q.x/y() - 2q.z() < 0 ? then z > x/y
                    if (projection(comp, comp) - projection(2, 2) > 0.0) {
                        comp = 2;
                    }

                    if (comp == 1)
                        tilt_color = Color3::red();
                    else if (comp == 2)
                        tilt_color = Color3::green();
                }

                // Apply the color
                for (int d = 0; d < 3; d++) {

                    verts[i].color[d] = tilt_color[d];

                    if (_vortexShell->invertNormals)
                        verts[i].normal[d] *= -1;
                }

            }

            LC_INFO("Successfully generated surface (verts = {0}, indices = {1})", nVert, nInd);

            // Fill magnum class with generated surface data
            _vortexShell->surface.Init((LC::Surface::Vertex*)verts, nVert, indices, nInd, _widget.preimage_translate);

            // Delete data
            delete[] verts;
            delete[] indices;

            _vortexShell->mesh = _vortexShell->surface.Mesh();

            // Add mesh to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_preimageManipulator, _phongShader, *(_vortexShell->mesh), _vortexShell->draw, _vortexShell->noCull, _transparentNormalDrawables };
        }
    }

    if (_baryonIsosurface && _baryon_density) {
        std::array<int, 3> voxBaryon = data->voxels;

        unsigned int baryVol = voxBaryon[0] * voxBaryon[1] * voxBaryon[2];
        
        color2[3] = _baryonIsosurface->alpha;
        std::unique_ptr<float[]> iso = std::unique_ptr<float[]>(new float[baryVol]);

        for (int i = 0; i < baryVol; i++) {
            // Diffmag for helical axis
            iso[i] = abs(_baryon_density[i] - _baryonIsosurface->query);
        }


        // Count the number of points within the preimage
        unsigned int domainPoints = 0;
        for (unsigned int i = 0; i < baryVol; i++) if (iso[i] < _baryonIsosurface->isoLevel) ++domainPoints;

        // If points found is greater than half of volume, then invert domain
        if (domainPoints > baryVol / 2) {
            LC_INFO("Preimage points found [{0}/{1}] exceeds half of volume: Inverting domain", domainPoints, baryVol);
            for (unsigned int i = 0; i < baryVol; i++) {
                if (iso[i] > _baryonIsosurface->isoLevel) iso[i] = 0.0f;
                else iso[i] = 10.0f;
            }
        }

        gen.GenerateSurface(iso.get(), _baryonIsosurface->isoLevel, voxBaryon, cell, color2);

        if (gen.isSurfaceValid()) {

            unsigned int nVert = gen.NumSurfaceVertices();
            unsigned int nInd = gen.NumSurfaceIndices();

            LC::Math::IsoVertex* verts = gen.ReleaseSurfaceVertices();
            unsigned int* indices = gen.ReleaseSurfaceIndices();

            // Go through vertices and color according to max magnitude pion componenet
            // Note needs to include preimage_translate since embedded volumes can be chosen!

            for (int i = 0; i < nVert; i++) {
                // extract indices from position
                int xx = (verts[i].position[0] / data->cell_dims[0] + 0.5) * (voxBaryon[0] - 1);
                int yy = (verts[i].position[1] / data->cell_dims[1] + 0.5) * (voxBaryon[1] - 1);
                int zz = (verts[i].position[2] / data->cell_dims[2] + 0.5) * (voxBaryon[2] - 1);

                // Get chi
                unsigned int idx = xx + yy * voxBaryon[0] + zz * voxBaryon[0] * voxBaryon[1];

                //std::cout << xx << " " << yy << " " << zz << std::endl;
                auto q = _quaternion_field[idx];


                LC::scalar max_comp = abs(q.x());
                Color3 tilt_color = Color3::blue();

                if (max_comp < abs(q.y())) {
                    tilt_color = Color3::red();
                    max_comp = abs(q.y());
                }
                if (max_comp < abs(q.z())) {
                    tilt_color = Color3::green();
                }

                // Apply the color
                for (int d = 0; d < 3; d++)
                    verts[i].color[d] = tilt_color[d];

            }

            LC_INFO("Successfully generated surface (verts = {0}, indices = {1})", nVert, nInd);

            // Fill magnum class with generated surface data
            _baryonIsosurface->surface.Init((LC::Surface::Vertex*)verts, nVert, indices, nInd, _widget.preimage_translate);

            // Delete data
            delete[] verts;
            delete[] indices;

            _baryonIsosurface->mesh = _baryonIsosurface->surface.Mesh();

            // Add mesh to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_preimageManipulator, _phongShader,
                *(_baryonIsosurface->mesh), _baryonIsosurface->draw,
                _baryonIsosurface->noCull, _transparentNormalDrawables };
        }
    }

    if (_quaternion_field && _pionPreimage) {
        std::array<int, 3> voxBaryon = data->voxels;
        unsigned int baryVol = voxBaryon[0] * voxBaryon[1] * voxBaryon[2];

        color2[3] = _pionPreimage->alpha;
        std::unique_ptr<float[]> iso = std::unique_ptr<float[]>(new float[baryVol]);

        Eigen::Quaterniond pi(0., _pionPreimage->pi1, _pionPreimage->pi2, _pionPreimage->pi3);

        for (int i = 0; i < baryVol; i++) {
            Eigen::Quaterniond q = _quaternion_field[i];

            // Get sin(theta/2)
            double pix = abs(q.x());
            double piy = abs(q.y());
            double piz = abs(q.z());
   
            iso[i] = q.w() * q.w() + pow(pix - abs(pi.x()), 2)
                + pow(piy - abs(pi.y()), 2)
                + pow(piz - abs(pi.z()), 2);
        }


        // Count the number of points within the preimage
        unsigned int domainPoints = 0;
        for (unsigned int i = 0; i < baryVol; i++) if (iso[i] < _pionPreimage->isoLevel) ++domainPoints;

        // If points found is greater than half of volume, then invert domain
        if (domainPoints > baryVol / 2) {
            LC_INFO("Preimage points found [{0}/{1}] exceeds half of volume: Inverting domain", domainPoints, baryVol);
            for (unsigned int i = 0; i < baryVol; i++) {
                if (iso[i] > _pionPreimage->isoLevel) iso[i] = 0.0f;
                else iso[i] = 10.0f;
            }
        }

        gen.GenerateSurface(iso.get(), _pionPreimage->isoLevel, voxBaryon, cell, color2);

        if (gen.isSurfaceValid()) {

            unsigned int nVert = gen.NumSurfaceVertices();
            unsigned int nInd = gen.NumSurfaceIndices();

            LC::Math::IsoVertex* verts = gen.ReleaseSurfaceVertices();
            unsigned int* indices = gen.ReleaseSurfaceIndices();

            // Go through vertices and color according to max magnitude pion componenet
            // Note needs to include preimage_translate since embedded volumes can be chosen!

            for (int i = 0; i < nVert; i++) {
                // extract indices from position
                int xx = (verts[i].position[0] / data->cell_dims[0] + 0.5) * (voxBaryon[0] - 1);
                int yy = (verts[i].position[1] / data->cell_dims[1] + 0.5) * (voxBaryon[1] - 1);
                int zz = (verts[i].position[2] / data->cell_dims[2] + 0.5) * (voxBaryon[2] - 1);

                // Get chi
                unsigned int idx = xx + yy * voxBaryon[0] + zz * voxBaryon[0] * voxBaryon[1];

                //std::cout << xx << " " << yy << " " << zz << std::endl;
                auto q = _quaternion_field[idx];


                LC::scalar max_comp = abs(q.x());
                Color3 tilt_color = Color3::blue();

                if (max_comp < abs(q.y())) {
                    tilt_color = Color3::red();
                    max_comp = abs(q.y());
                }
                if (max_comp < abs(q.z())) {
                    tilt_color = Color3::green();
                }

                // Apply the color
                for (int d = 0; d < 3; d++)
                    verts[i].color[d] = tilt_color[d];

            }

            LC_INFO("Successfully generated surface (verts = {0}, indices = {1})", nVert, nInd);

            // Fill magnum class with generated surface data
            _pionPreimage->surface.Init((LC::Surface::Vertex*)verts, nVert, indices, nInd, _widget.preimage_translate);

            // Delete data
            delete[] verts;
            delete[] indices;

            _pionPreimage->mesh = _pionPreimage->surface.Mesh();

            // Add mesh to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_preimageManipulator, _phongShader,
                *(_pionPreimage->mesh), _pionPreimage->draw,
                _pionPreimage->noCull, _transparentNormalDrawables };
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