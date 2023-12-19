#include "Widget.h"
#include "find_components.h"
#include "ZProfile.h"
#include "LehmanCluster.h"
#include "CrossSection.h"
#include "SmoothIsosurface.h"
#include "QuaternionFromBasis.h"
#include <tuple>
#include <filesystem>
#include <set>

using FOFDSolver = LC::FrankOseen::Electric::FOFDSolver;
using Dataset = FOFDSolver::Dataset;

#define MAX_GRAPH_POINTS 250
#define EXTENDEDMC 1
// Current use case for interpolating preimages is inefficient,
// use to interpolate director field THEN copy before passing to isoGenerator
#define TESTINTERP 0

std::tuple<std::vector<Eigen::Vector3d>, std::vector<int>> exportVortexKnot(const std::vector<MeshLib::PNCVertex<float>>& vertices,
    const std::vector<MeshLib::Triangle>& triangles,
    LC::scalar min_point_dist, const std::string& fname);

std::vector<Eigen::Vector3d> detect_heliknoton_COM(std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>>&,const std::vector<Eigen::Vector3d>&);

struct BoxInstanceData {
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

struct HopfDensityIsosurface {

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

struct NematicMultiplaneManager {
    struct NematicPlane {
        NematicPlane(int id = -1) : node_density(35),
                                    sheet_density(100),
                                    plane_normal(Magnum::Vector2{ M_PI/2.f, 3.f * M_PI / 2.f }),
                                    plane_id(id),
                                    plane_diam(0.222f),
                                    plane_center(Magnum::Vector3{0.f, 0.f, 0.f}),
                                    draw(false),
                                    opentab(true),
                                    plane_preview(true),
                                    initialized(false),
                                    vortexTolerance(0.0f),
                                    scaling(0.193f),
                                    rod_length(3.3f),
                                    offset(-0.02f),
                                    pos_gamma(15.f),
                                    neg_gamma(15.f),
                                    cutoff(.2f),
                                    field(0) {
        }
        LC::NematicArray plane;
        std::unique_ptr<LC::DynamicColorSheet> sheet;
        GL::Mesh sheet_mesh;
        Magnum::Vector2 plane_normal; // (theta, phi) coordinates
        Magnum::Vector3 plane_center;
        int node_density;
        int sheet_density;
        float plane_diam;
        bool initialized;
        bool draw;
        bool opentab;
        int plane_id;
        int field; // 0 -> lambda, 1 -> chi, 2 -> tau
        Float offset; // Plane offset from rods
        bool plane_preview;
        Float scaling;
        Float rod_length;
        Float pos_gamma;
        Float neg_gamma;
        Float cutoff;
        float vortexTolerance;
    };

    NematicMultiplaneManager() {
        planes.reserve(10);
    }

    int Dim() const { return planes.size(); }

    void CleanAll() {
        planes.clear();
        selected_plane = -1;
    }

    void Draw(const Matrix4 &viewMatrix, const Matrix4 &projection, Widget &widget) {
        // Iterate through the planes and draw
        for (auto& plane : planes) {

            plane.plane.specular = widget.specular;
            plane.plane.diffuse = widget.diffuse;
            plane.plane.ambient = widget.ambient;
            plane.plane.selected_drawType = LC::NematicArray::DrawType::Cylinder;

            if (plane.draw && plane.initialized)
                plane.plane.Draw(viewMatrix,projection);
        }
    }

    void Update(Dataset *data) {
        // Update planes
        LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
        LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
        LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);
        uint Nvox = data->voxels[0] * data->voxels[1] * data->voxels[2];

        for (auto& plane : planes) {
            
            if (plane.initialized) {
                for (int i = 0; i < plane.plane.polyPositions.size(); i++) {
                    auto position = plane.plane.polyPositions[i];

                    Float scale = plane.plane.scale / (Float)pow(plane.plane.polyPositions.size(), 1.0f / 3.0f);

                    // Translate position coords to voxel space
                    int xi = (data->cell_dims[0] * 0.5 + position[0]) / dx;
                    int yi = (data->cell_dims[1] * 0.5 + position[1]) / dy;
                    int zi = (data->cell_dims[2] * 0.5 + position[2]) / dz;

                    // Make sure within bounds
                    if (xi < 0 || yi < 0 || zi < 0 ||
                        xi >= data->voxels[0] || yi >= data->voxels[1] || zi >= data->voxels[2])
                        continue;

                    uint id = xi + data->voxels[0] * yi + data->voxels[1] * data->voxels[2] * zi;


                    // Get theta and phi for director
                    LC::scalar fx, fy, fz;
                    // Field is set to director
                    if (plane.field == 0) {
                        fx = data->directors[id];
                        fy = data->directors[id + Nvox];
                        fz = data->directors[id + 2 * Nvox];
                    }
                    else if (plane.field == 1 || plane.field == 2) {
                        // Field is set to chi
                        // Compute the chi field at this point
                        // If the point is too close to the bulk, just set to homeotropic
                        auto chi_field = LC::Math::ChiralityField_local(data->directors.get(), xi, yi, zi, data->voxels, data->cell_dims, (LC::scalar)plane.vortexTolerance);
                        fx = chi_field(0);
                        fy = chi_field(1);
                        fz = chi_field(2);
                    }

                    float theta = acos(fz);
                    float phi = atan2(fy, fx);

                    // Will depend on the chosen field to draw

                    Magnum::Color3 director_color;
                    if (plane.field == 0) {
                        director_color = LC::Imaging::Colors::RungeSphere(theta, phi);
                    }
                    else if (plane.field == 1) {
                        float thetap = theta >= 0.5f * M_PI ? M_PI - theta : theta;
                        director_color = LC::Imaging::Colors::RungeSphere(thetap, 2 * phi);
                    }
                    else if (plane.field == 2) {
                        // Black directors
                        // and splay bend colors;
                        director_color = Color3(rod_color[0], rod_color[1], rod_color[2]);

                    }

                    // Apply rotation based on director
                    float hPi = M_PI / 2.f;
                    Matrix4 transformation = Matrix4::translation(position)
                        * Matrix4::rotationZ(Rad{ -hPi + float(phi) })
                        * Matrix4::rotationX(Rad{ hPi - float(theta) })
                        * Matrix4::scaling({ scale,scale,scale });

                    plane.plane.polyInstanceData[i].transformationMatrix = transformation;
                    plane.plane.polyInstanceData[i].color = director_color;
                    plane.plane.polyInstanceData[i].normalMatrix = transformation.normalMatrix();
                }

                // Record values so that analysis can be done
                std::vector<Float> Ssb_plane(plane.sheet->vertices.size());
                Float Smin(0.f), Smax(0.f), Spos_avg(0.f), Sneg_avg(0.f);
                int NSpos = 0;
                int NSneg = 0;
                // Create a lambda func to access a particular node
                auto chi_f = [&](int xi, int yi, int zi) {
                    return LC::Math::ChiralityField_local(data->directors.get(), xi, yi, zi, data->voxels, data->cell_dims, (LC::scalar)plane.vortexTolerance);
                };

                auto qa_f = [](const Eigen::Vector3d& chi) {
                    return 3. * chi(0) * chi(0) - 1.;
                };
                auto qb_f = [](const Eigen::Vector3d& chi) {
                    return 3. * chi(0) * chi(1);
                };
                auto qc_f = [](const Eigen::Vector3d& chi) {
                    return 3. * chi(0) * chi(2);
                };
                auto qd_f = [](const Eigen::Vector3d& chi) {
                    return 3. * chi(1) * chi(1) - 1.;
                };
                auto qe_f = [](const Eigen::Vector3d& chi) {
                    return 3. * chi(1) * chi(2);
                };

                auto splay_bend = [&](int xi, int yi, int zi) {
                    auto _000 = chi_f(xi, yi, zi);
                    auto _p00 = chi_f(xi + 1, yi, zi);
                    auto _m00 = chi_f(xi - 1, yi, zi);
                    auto _0p0 = chi_f(xi, yi + 1, zi);
                    auto _0m0 = chi_f(xi, yi - 1, zi);
                    auto _00p = chi_f(xi, yi, zi + 1);
                    auto _00m = chi_f(xi, yi, zi - 1);
                    auto _pp0 = chi_f(xi + 1, yi + 1, zi);
                    auto _pm0 = chi_f(xi + 1, yi - 1, zi);
                    auto _mp0 = chi_f(xi - 1, yi + 1, zi);
                    auto _mm0 = chi_f(xi - 1, yi - 1, zi);
                    auto _p0p = chi_f(xi + 1, yi, zi + 1);
                    auto _m0p = chi_f(xi - 1, yi, zi + 1);
                    auto _p0m = chi_f(xi + 1, yi, zi - 1);
                    auto _m0m = chi_f(xi - 1, yi, zi - 1);
                    auto _0pp = chi_f(xi, yi + 1, zi + 1);
                    auto _0pm = chi_f(xi, yi + 1, zi - 1);
                    auto _0mp = chi_f(xi, yi - 1, zi + 1);
                    auto _0mm = chi_f(xi, yi - 1, zi - 1);

                    Float qa200 = (qa_f(_p00) + qa_f(_m00) - qa_f(_000)) / (dx * dx);
                    Float qa002 = (qa_f(_00p) + qa_f(_00m) - qa_f(_000)) / (dz * dz);
                    Float qd002 = (qd_f(_00p) + qd_f(_00m) - qd_f(_000)) / (dz * dz);
                    Float qd020 = (qd_f(_0p0) + qd_f(_0m0) - qd_f(_000)) / (dy * dy);
                    Float qb110 = (qb_f(_pp0) + qb_f(_mm0) - qb_f(_pm0) - qb_f(_mp0)) / (4.f * dx * dy);
                    Float qc101 = (qc_f(_p0p) + qc_f(_m0m) - qc_f(_p0m) - qc_f(_m0p)) / (4.f * dx * dz);
                    Float qe011 = (qe_f(_0pp) + qe_f(_0mm) - qe_f(_0pm) - qe_f(_0mp)) / (4.f * dy * dz);

                    return -qa002 + qa200 + 2.f * (qb110 + qc101) - qd002 + qd020 + 2.f * qe011;
                };

                for (int i = 0; i < plane.sheet->vertices.size(); i++) {
                    auto position = plane.sheet->vertices[i].position;

                    // Translate position coords to voxel space
                    int xi = (data->cell_dims[0] * 0.5 + position[0]) / dx;
                    int yi = (data->cell_dims[1] * 0.5 + position[1]) / dy;
                    int zi = (data->cell_dims[2] * 0.5 + position[2]) / dz;

                    // Make sure within bounds
                    if (xi < 0 || yi < 0 || zi < 0 ||
                        xi >= data->voxels[0] || yi >= data->voxels[1] || zi >= data->voxels[2])
                        continue;

                    uint id = xi + data->voxels[0] * yi + data->voxels[1] * data->voxels[2] * zi;


                    // Get theta and phi for director
                    LC::scalar fx(1.f), fy(0.f), fz(0.f);
                    // Field is set to director
                    if (plane.field == 0) {
                        fx = data->directors[id];
                        fy = data->directors[id + Nvox];
                        fz = data->directors[id + 2 * Nvox];
                    }
                    else if (plane.field == 1) {
                        // Field is set to chi
                        // Compute the chi field at this point
                        // If the point is too close to the bulk, just set to homeotropic
                        auto chi_field = chi_f(xi, yi, zi);

                        fx = chi_field(0);
                        fy = chi_field(1);
                        fz = chi_field(2);
                    }
                    else if (plane.field == 2) {

                        Float cijk[2][2][2];
                        Float cjk[2][2];
                        Float ck[2];

                        // Use trilinear interpolation
                        Float x = (data->cell_dims[0] * 0.5 + position[0]);
                        Float y = (data->cell_dims[1] * 0.5 + position[1]);
                        Float z = (data->cell_dims[2] * 0.5 + position[2]);
                        // Bottom corner points
                        int i_0 = (data->cell_dims[0] * 0.5 + position[0]) / dx;
                        int j_0 = (data->cell_dims[1] * 0.5 + position[1]) / dy;
                        int k_0 = (data->cell_dims[2] * 0.5 + position[2]) / dz;
                        Float x_0 = i_0 * dx;
                        Float y_0 = j_0 * dy;
                        Float z_0 = k_0 * dz;
                        // Interpolation coordinates
                        Float x_d = (x - x_0) / dx;
                        Float y_d = (y - y_0) / dy;
                        Float z_d = (z - z_0) / dz;
                        // Fill the interpolant
                        for (int a = 0; a < 2; a++)
                            for (int b = 0; b < 2; b++)
                                for (int c = 0; c < 2; c++)
                                    cijk[a][b][c] = splay_bend(i_0 + a, j_0 + b, k_0 + c);

                        // Begin interpolation
                        for (int b = 0; b < 2; b++)
                            for (int c = 0; c < 2; c++)
                                cjk[b][c] = (1.f - x_d) * cijk[0][b][c] + x_d * cijk[1][b][c];

                        // Begin interpolation
                        for (int c = 0; c < 2; c++)
                            ck[c] = (1.f - y_d) * cjk[0][c] + y_d * cjk[1][c];


                        //Ssb_plane[i] = splay_bend(xi, yi, zi);
                        Ssb_plane[i] = (1.f - z_d) * ck[0] + z_d * ck[1];

                        if (Smin > Ssb_plane[i])
                            Smin = Ssb_plane[i];
                        if (Smax < Ssb_plane[i])
                            Smax = Ssb_plane[i];

                        if (Ssb_plane[i] > 0.f) {
                            Spos_avg += Ssb_plane[i];
                            NSpos++;
                        }
                        if (Ssb_plane[i] < 0.f) {
                            Sneg_avg += Ssb_plane[i];
                            NSneg++;
                        }

                    }

                    float theta = acos(fz);
                    float phi = atan2(fy, fx);

                    Magnum::Color3 director_color;
                    if (plane.field == 0) {
                        director_color = LC::Imaging::Colors::RungeSphere(theta, phi);
                        plane.sheet->vertices[i].color = { director_color, 0.9f };
                    }
                    else if (plane.field == 1) {
                        float thetap = theta >= 0.5f * M_PI ? M_PI - theta : theta;
                        director_color = LC::Imaging::Colors::RungeSphere(thetap, 2 * phi);
                        plane.sheet->vertices[i].color = { director_color, 0.9f };
                    }
                    
                }
                if (plane.field == 2) {

                    if (NSpos > 0)
                        Spos_avg /= NSpos;
                    if (NSneg > 0)
                        Sneg_avg /= NSneg;

                    Magnum::Color3 blue(0.f, 0.f, 1.f);
                    Magnum::Color3 neutral(0.0f, 0.0f, 0.0f);
                    Magnum::Color3 yellow(1.f, 1.f, 0.f);
                    for (int i = 0; i < plane.sheet->vertices.size(); i++) {

                        // Set the color
                        Magnum::Color3 director_color;
                        if (Ssb_plane[i] < 0.f) {
                            // Parametrize the color from 0 to the most negative splay bend
                            Float t = 0.f; Float t0 = 0.5f; Float tc = 0.75f;
                            if (Smin < 0.f) {
                                t = Ssb_plane[i] / Smin;
                                t0 = Sneg_avg / Smin;
                                // Choose cutoff between t0 and 1
                                Float u = plane.cutoff;
                                tc = (1.f - u) * t0 + u;
                            }
                            t = 1. / (1. + exp(-(t - tc) * plane.neg_gamma));
                            director_color = (1.f - tc) * neutral + t * yellow;
                        }
                        else {
                            // Parametrize the color from 0 to the most positive splay bend
                            Float t = 0.f; Float t0 = 0.5f; Float tc = 0.75f;
                            if (Smax > 0.f) {
                                t = Ssb_plane[i] / Smax;
                                // Normalize avg
                                t0 = Spos_avg / Smax;
                                // Choose cutoff between t0 and 1
                                Float u = plane.cutoff;
                                tc = (1.f - u) * t0 + u;
                            }
                            t = 1. / (1. + exp(-(t - tc) * plane.pos_gamma));
                            director_color = (1.f - t) * neutral + t * blue;
                        } 
                        plane.sheet->vertices[i].color = { director_color, 0.9f };
                    }

                    // Apply a Gaussian blending
                    std::vector<Color4> colors(plane.sheet->vertices.size());
                    // Copy the colors above
                    for (int i = 0; i < plane.sheet->vertices.size(); i++) {
                        colors[i] = plane.sheet->vertices[i].color;
                    }
                    // Box blur the colors
                    for (int i = 1; i < plane.sheet->NX-1; i++) {
                        for (int j = 1; j < plane.sheet->NY-1; j++) {
                            Color4 sum(0.f, 0.f, 0.f, 0.f);
                            for (int x = -1; x <= 1; x++)
                                for (int y = -1; y <= 1; y++)
                                    sum = sum + colors[i + x + plane.sheet->NX * (j + y)];
                            plane.sheet->vertices[i + plane.sheet->NX * j].color = sum / 9.f;
                        }
                    }

                }

                // Update sheet
                Trade::MeshData meshData = plane.sheet->Data();
                plane.sheet->vertexBuffer.setData(MeshTools::interleave(meshData.positions3DAsArray(), meshData.colorsAsArray()), GL::BufferUsage::DynamicDraw);
            }
        }
    }

    void Gui(Dataset *data,
        LC::Drawable::Object3D *manipulator,
        Shaders::VertexColor3D &transparentShader,
        SceneGraph::DrawableGroup3D &transparentDrawables) {

        if (ImGui::Button("Add new plane")) {
            NewPlane();
        }

        int current_plane_size;

        do {
            current_plane_size = planes.size();
            for (int i = 0; i < current_plane_size; i++) {
                if (!planes[i].opentab) {
                    // Stop drawing plane and remove
                    planes[i].draw = false;
                    planes.erase(planes.begin() + i);
                    break;
                }
            }

        } while (current_plane_size != planes.size());

        ImGui::ColorEdit4("Rod color##Multiplane-viewer", &rod_color[0], ImGuiColorEditFlags_NoInputs);

        if (planes.size()) {
            ImGui::Text("Planes");
            ImGui::Separator();
        }
        char txtbuffer[30];
        float maxdim = 0.5f * (std::max)({ cell[0], cell[1], cell[2] });

        if (ImGui::BeginTabBar("Selected planes", ImGuiTabBarFlags_Reorderable | ImGuiTabBarFlags_AutoSelectNewTabs | ImGuiTabBarFlags_FittingPolicyScroll)) {
            int ctr = 0;
            for (auto& p : planes) {

                int n2 = sprintf(txtbuffer, "Plane %d", ++ctr);
                // Copy to string
                std::string txt;
                for (int i = 0; i < n2; i++) txt += txtbuffer[i];

                if (p.opentab && ImGui::BeginTabItem(txt.c_str(), &p.opentab, ImGuiTabItemFlags_None))
                {
                    ImGui::Text("%s", txt.c_str());
                    ImGui::Separator();
                    // Plane parameters
                    ImGui::RadioButton(("lambda##plane" + txt).c_str(), &p.field, 0);
                    ImGui::SameLine();
                    ImGui::RadioButton(("chi##plane" + txt).c_str(), &p.field, 1);
                    ImGui::SameLine();
                    ImGui::RadioButton(("chi-SB##plane" + txt).c_str(), &p.field, 2);
                    if (p.field == 1 || p.field == 2) {
                        ImGui::SliderFloat(("Vortex tolerance##plane" + txt).c_str(), &p.vortexTolerance, 0.0f, 0.1f);
                    }
                    ImGui::SliderFloat(("SB (+) gamma##plane" + txt).c_str(), &p.pos_gamma, 0.0f, 15.f);
                    ImGui::SliderFloat(("SB (-) gamma##plane" + txt).c_str(), &p.neg_gamma, 0.0f, 15.f);
                    ImGui::SliderFloat(("SB cutoff##plane" + txt).c_str(), &p.cutoff, 0.0f, 1.0f);
                    ImGui::Text("Densities");
                    ImGui::SliderInt(("Density##plane" + txt).c_str(), &p.node_density, 5, 60);
                    ImGui::SliderInt(("Sheet density##plane" + txt).c_str(), &p.sheet_density, 20, 150);
                    ImGui::SliderFloat(("Scaling##plane" + txt).c_str(), &p.scaling, 0.1f, 2.f);
                    ImGui::SliderFloat(("Rod length##plane" + txt).c_str(), &p.rod_length, 2.f, 6.f);
                    
                    ImGui::Separator();
                    ImGui::Text("Plane Orientation");
                    ImGui::SliderFloat(("Theta##plane" + txt).c_str(), &p.plane_normal[0], 0.0f, M_PI);
                    ImGui::SliderFloat(("Phi##plane" + txt).c_str(), &p.plane_normal[1], 0.0f, 2.f*M_PI);
                    ImGui::SliderFloat(("Offset##plane" + txt).c_str(), &p.offset, -0.1f, 0.1f);
                    ImGui::Separator();
                    ImGui::Text("Plane position");
                    ImGui::SliderFloat(("x##plane" + txt).c_str(), &p.plane_center[0], -0.5f * cell[0], 0.5f * cell[0]);
                    ImGui::SliderFloat(("y##plane" + txt).c_str(), &p.plane_center[1], -0.5f * cell[1], 0.5f * cell[1]);
                    ImGui::SliderFloat(("z##plane" + txt).c_str(), &p.plane_center[2], -0.5f * cell[2], 0.5f * cell[2]);
                    ImGui::Separator();
                    ImGui::Text("Plane dimensions");
                    ImGui::SliderFloat(("D##plane" + txt).c_str(), &p.plane_diam, 0, maxdim);
                    ImGui::Separator();
                    ImGui::Checkbox(("Draw##plane" + std::to_string(p.plane_id)).c_str(), &p.draw);
                    ImGui::SameLine();
                    ImGui::Checkbox("Auto-orient", &p.plane_preview);
                    ImGui::SameLine();
                    if (ImGui::Button("Init")) {

                        // Set initialized variable
                        p.initialized = true;

                        // Initialize the plane!
                        std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* data, std::size_t i) {
                            Magnum::Vector3* pos = (Magnum::Vector3*)data;
                            return pos[i];
                        };
                        
                        int density = p.node_density;
                        uint Nr = p.plane_diam * density;
                        uint Nr_s = p.plane_diam * p.sheet_density;
                        Float dr = p.plane_diam / (Nr - 1);
                        Float dr_s = p.plane_diam / (Nr_s - 1);
                        uint sz = 0;

                        // Sample points on a plane
                        std::vector<Vector2> sampled_points;
                        for (int i = 0; i < Nr; i++) {
                            for (int j = 0; j < Nr; j++) {

                                // Sampling with spiral
                                //Float ii = (i + j * Nr) / (float)Nr / (float)Nr;
                                //Float u = p.plane_diam * sqrt(ii) * cos(1.618*(i + j * Nr) * 2.f * M_PI);
                                //Float v = p.plane_diam * sqrt(ii) * sin(1.618*(i + j * Nr) * 2.f * M_PI);
                                // Regular sampling
                                Float u = i * dr - p.plane_diam * 0.5f;
                                Float v = j * dr - p.plane_diam * 0.5f;
                                sampled_points.emplace_back(u, v);
                            }
                        }
                        sz = sampled_points.size();

                        // Initialize position data
                        std::unique_ptr<Vector3[]> position_data(new Vector3[sz]);
                        std::unique_ptr<Vector3[]> sheet_position_data(new Vector3[Nr_s * Nr_s]);

                        // Specify the normal vector
                        Float sth = sin(p.plane_normal[0]);
                        Float cth = cos(p.plane_normal[0]);
                        Float sph = sin(p.plane_normal[1]);
                        Float cph = cos(p.plane_normal[1]);
                        Eigen::Vector3f nhat(sth * cph, sth * sph, cth);

                        // Build an orthogonal set of vectors to nhat
                        Eigen::Vector3f e1(1.f, 0.f, 0.f), e2;

                        if(abs(e1.dot(e2)) > 0.99f) {
                            e1(1) = 1.f;
                        }
                        // Orthogonalize e1 wrt nhat
                        e1 = e1.cross(nhat);
                        e1.normalize();
                        e2 = nhat.cross(e1);

                        for (uint i = 0; i < sz; i++) {
                            auto point = sampled_points[i];
                            Float u = point[0];
                            Float v = point[1];

                            position_data[i][0] = u * e1(0) + v * e2(0) + p.plane_center[0];
                            position_data[i][1] = u * e1(1) + v * e2(1) + p.plane_center[1];
                            position_data[i][2] = u * e1(2) + v * e2(2) + p.plane_center[2];
                        }

                        for (uint i = 0; i < Nr_s; i++) {
                            for (uint j = 0; j < Nr_s; j++) {
                                Float u = 1.2f * (dr_s * i - p.plane_diam * 0.5f);
                                Float v = 1.2f * (dr_s * j - p.plane_diam * 0.5f);
                                sheet_position_data[i + j * Nr_s][0] = u * e1(0) + v * e2(0) + p.plane_center[0] + p.offset * nhat(0);
                                sheet_position_data[i + j * Nr_s][1] = u * e1(1) + v * e2(1) + p.plane_center[1] + p.offset * nhat(1);
                                sheet_position_data[i + j * Nr_s][2] = u * e1(2) + v * e2(2) + p.plane_center[2] + p.offset * nhat(2);
                            }
                        }

                        p.plane.scale = p.scaling/sqrt(density)/1.618;
                        p.plane.hLength = p.rod_length;
                        p.plane.Init((void*)position_data.get(), access, sz);

                        //Initialize the sheet
                        p.sheet = std::make_unique<LC::DynamicColorSheet>();
                        p.sheet->Set(Nr_s, Nr_s, p.plane_diam, p.plane_diam);


                        std::function<Vector3(void*,std::size_t)> sfunc = [&](void*dat,std::size_t i) {
                            Vector3* dat_vec = static_cast<Vector3*>(dat);
                            return dat_vec[i];
                        };

                        p.sheet->InitManual((void*)sheet_position_data.get(), sfunc);
                        p.sheet_mesh = p.sheet->Mesh();
                        new LC::Drawable::TransparentFlatDrawable{ *manipulator, transparentShader,
                            p.sheet_mesh, p.draw, Vector3{0.0f, 0.0f, 0.0f}, transparentDrawables };
                        
                        
                        // Call update
                        Update(data);
                    }

                    ImGui::EndTabItem();
                }
            }
            ImGui::EndTabBar();
        }
    }
    void NewPlane() {
        // Add a new plane
        planes.push_back(NematicPlane(global_plane_counter++));
        if (planes.size() == 1)
            selected_plane = 0;
    }

    std::vector<NematicPlane> planes;
    std::array<LC::scalar, 3> cell;
    Color4 rod_color = { 0.5f, 0.5, 0.5f, 1.f };
    int selected_plane = -1;
    
    int global_plane_counter = 0; // Count how many planes have been created
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
    void drawBoxes();
    void imageMenu();
    void handlePOMWindow();
    void handlePreimageWindow();
    void handleVortexKnotWindow();
    void handleKnotInteractionWindow();
    void nematicPreimageManager();
    void handleNonlinearImagingWindow();
    void handleLCINFOWindow();
    void handleModificationWindow();
    void handleZProfileWindow();
    void generatePionTriplet();
    void handleMultiplaneWindow();
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
    void findVortexKnotComponents2(bool saveObj = false, const std::string& obj_name = "test.obj");

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
    std::list<S2Fiber> _nematicPreimages;
    std::list<VortexKnot> _vortexKnot, _processedVortexKnot;
    Containers::Optional<VortexShell> _vortexShell;
    Containers::Optional<PionPreimage> _pionPreimage;
    Containers::Optional<BaryonIsosurface> _baryonIsosurface;
    Containers::Optional<HopfDensityIsosurface> _hopfionDensityIsosurface;
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
    (*data).Configuration(Dataset::Heliknoton(Q,data->voxels,data->cell_dims,1.,1.));

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

        Dataset::Config helis = Dataset::Heliknoton(1,data->voxels,data->cell_dims, translations, 0.50, 3);

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
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);


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

                    // Normalize directors
                    size_t grid_vol = dataloc->voxels[0] * dataloc->voxels[1] * dataloc->voxels[2];
                    size_t invalid_dirs = 0;
                    for (size_t i = 0; i < grid_vol; i++) {

                        LC::scalar mag = sqrt(dataloc->directors[i] * dataloc->directors[i] + dataloc->directors[i + grid_vol] * dataloc->directors[i + grid_vol]
                            + dataloc->directors[i + 2 * grid_vol] * dataloc->directors[i + 2 * grid_vol]);

                        if (mag != 0.)
                            for (int d = 0; d < 3; d++)
                                dataloc->directors[i + d * grid_vol] /= mag;
                        else
                            invalid_dirs++;
                    }

                    if (invalid_dirs)
                        LC_WARN("Loaded degenerate director field");

                    _widget.energy_series = std::list<LC::precision_scalar>{};
                    _widget.series_x_axis = std::list<LC::precision_scalar>{};
                    for (int i = 0; i < _widget.energy_series_vec.size(); i++) {
                        _widget.energy_series_vec[i] = 0.0;
                        _widget.series_x_axis_vec[i] = dataloc->numIterations + i;
                    }

                    // ** Only intended to be used for POM recording
                    _image_series = LC::Imaging::ImageSeries((std::int32_t)dataloc->voxels[0], (std::int32_t)dataloc->voxels[1], _header.readFile + "-POM");
                    for (int d = 0; d < 3; d++) _widget.celldims[d] = dataloc->cell_dims[d];

                };
                saveMenu(_widget.updateImageFromLoad, loadAction);
            }

            // Background color
            ImGui::ColorEdit4("Background Color##BackgroundColor", &_widget.clearColor[0], ImGuiColorEditFlags_NoInputs);

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
                handleVortexKnotWindow();
                handleKnotInteractionWindow();
                handleMultiplaneWindow();

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
                ImGui::SameLine();
                ImGui::Checkbox("S2 colors", &_widget.S2colors);


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
            static bool s_flip_manipulator_rot = false;
            Matrix4 inversion_mat = Matrix4::rotationX(Rad(s_flip_manipulator_rot * M_PI));

            ImGui::Text("Go to");
            ImGui::SameLine();
            if (ImGui::Button("xy")) {
                _manipulator->setTransformation(inversion_mat);
            }
            ImGui::SameLine();
            if (ImGui::Button("xz")) {
                _manipulator->setTransformation(inversion_mat * Matrix4::rotationZ(Rad(M_PI)) * Matrix4::rotationY(Rad(M_PI)) * Matrix4::rotationX(Rad(M_PI / 2.0f)));
            }
            ImGui::SameLine();
            if (ImGui::Button("yz")) {
                _manipulator->setTransformation(inversion_mat * Matrix4::rotationZ(Rad(-M_PI / 2.0f)) * Matrix4::rotationY(Rad(-M_PI / 2.0f)));
            }
            ImGui::SameLine();
            ImGui::Checkbox("Flip##Manipulator", &s_flip_manipulator_rot);

            ImGui::PushItemWidth(75.f);
            ImGui::InputFloat("Specular", &_widget.specular);
            ImGui::SameLine();
            ImGui::InputFloat("Diffuse", &_widget.diffuse);
            ImGui::SameLine();
            ImGui::InputFloat("Ambient", &_widget.ambient);
            ImGui::PopItemWidth();

            
            dropDownMenu<LC::NematicArray::DrawType>("Nematic draw type", _widget.nematicDrawType, LC::NematicArray::DrawMap());
            

            if (ImGui::Button("Update Image")  || _widget.updateImageFromLoad )
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

            // Compute change in free en:
            LC::scalar difference = 2.0 * _widget.energyErrorThreshold;
            if (_widget.energy_series.size() > 1) {
                std::list<LC::precision_scalar>::iterator it_end_m1 = _widget.energy_series.begin();
                std::advance(it_end_m1, _widget.energy_series.size() - 2);
                std::list<LC::precision_scalar>::iterator it_end = it_end_m1;
                std::advance(it_end, 1);
                
                difference = *it_end - *it_end_m1;
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
        if (_crossSections[id].draw_nematic && _crossSections[id].nematic) {
            _crossSections[id].nematic->specular = _widget.specular;
            _crossSections[id].nematic->diffuse = _widget.diffuse;
            _crossSections[id].nematic->ambient = _widget.ambient;
            _crossSections[id].nematic->selected_drawType = _widget.nematicDrawType;
            _crossSections[id].nematic->Draw(_camera->cameraMatrix() * _manipulator->transformationMatrix(), _camera->projectionMatrix());
        }
    }

    // Draw planes from plane selector
    _multiplane_manager.Draw(_camera->cameraMatrix()* _manipulator->transformationMatrix(), _camera->projectionMatrix(), _widget);

    // Draw vortex points
    for (auto iter = _processedVortexKnot.begin(); iter != _processedVortexKnot.end(); iter++) {
        if (iter->draw)
            iter->points.Draw(_camera->cameraMatrix()* _manipulator->transformationMatrix(), _camera->projectionMatrix());
    }


    // Draw z profile
    if (_widget.drawProfile && _zprofile) {
        _zprofile->Draw(_camera->cameraMatrix() * _manipulator->transformationMatrix(), _camera->projectionMatrix());
    }

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
        _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, true);
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
        _relaxFuture.first = LC::Solver::RelaxAsync(_solver.get(), _widget.cycle, std::launch::async, true);
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

    LC::Math::ChiralityField(nn.get(), chi_field, data->voxels, { (float)data->cell_dims[0],(float)data->cell_dims[1],(float)data->cell_dims[2] }, valid_field, true, _widget.isoLevel);
    
    // Extract pion triplet from n and chi
    /*
        Reference triplet: lambda = (1, 0, 0)
                           tau    = (-1, 0, 0)
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
                else {
                    // Use the convention that chi_z must always be positive
                    if (chi_field[full_idx + 2 * vol] < 0.)
                        chi = Eigen::Vector3f{ -chi_field[full_idx], -chi_field[full_idx + vol], -chi_field[full_idx + 2 * vol] };
                    else
                        chi = Eigen::Vector3f{ chi_field[full_idx], chi_field[full_idx + vol], chi_field[full_idx + 2 * vol] };

                }

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
                if (!_widget.POM || ax != Axis::z) {
                    _crossSections[id].section.second->vertices[cidx].color = { director_color, alpha };
                }
                if (_widget.S2colors) _crossSections[id].nematic->polyInstanceData[cidx_red].color = director_color;
                else _crossSections[id].nematic->polyInstanceData[cidx_red].color = Color3{ 0.5f, 0.5f, 0.5f };

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
                //_crossSections[id].nematic->polyInstanceData[cidx_red].normalMatrix = transformation.normalMatrix();
            }
        }

        // Update sheet
        Trade::MeshData meshData = _crossSections[id].section.second->Data();
        _crossSections[id].section.second->vertexBuffer.setData(MeshTools::interleave(meshData.positions3DAsArray(), meshData.colorsAsArray()), GL::BufferUsage::DynamicDraw);



    }

}

void Sandbox::POM(std::string override_name) {
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
        std::string name;
        if (override_name.empty())
            name = getLoadFile();
        else
            name = override_name;

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
    _multiplane_manager.CleanAll();

    for (int id = 0; id < 3; id++) {

        Axis ax = static_cast<Axis>(id);

        int i = (ax == Axis::x) ? 1 : 0;
        int j = (ax == Axis::y) ? 2 : i + 1;

        // Start by only drawing the xy plane
        if (ax != Axis::z) {
            _crossSections[id].draw = false;
            _crossSections[id].draw_nematic = false;
        }

        int next = _crossSections[id].nematic_reduction_factor;
        _crossSections[id].nematic = LC::NematicArray{};
        unsigned int size = data->voxels[i] * data->voxels[j];
        int vox_red_x = (data->voxels[i] + next - 1) / next;
        int vox_red_y = (data->voxels[j] + next - 1) / next;
        unsigned int size_reduced = vox_red_x * vox_red_y;
        std::unique_ptr<Vector3[]> positions(new Vector3[size_reduced]);
        float dx = data->cell_dims[i] / (data->voxels[i] - 1);
        float dy = data->cell_dims[j] / (data->voxels[j] - 1);

        for (int x = 0; x < data->voxels[i]; x+=next) {
            for (int y = 0; y < data->voxels[j]; y+=next) {

                LC::scalar x_red_coord = x / (LC::scalar)(data->voxels[i]-1); // [0,1]
                LC::scalar y_red_coord = y / (LC::scalar)(data->voxels[j]-1); // [0,1]
                
                int x_red = ceil(x_red_coord * (vox_red_x - 1.1));
                int y_red = ceil(y_red_coord * (vox_red_y - 1.1));

                positions[x_red + y_red * vox_red_x] = { 0.f, 0.f, 0.f };
                positions[x_red + y_red * vox_red_x][i] = -data->cell_dims[i] * 0.5 + x * dx;
                positions[x_red + y_red * vox_red_x][j] = -data->cell_dims[j] * 0.5 + y * dy;
            }
        }
        // Now we can initialize the sphere array
        std::function<Magnum::Vector3(void*, std::size_t)> access = [](void* data, std::size_t i) {
            Magnum::Vector3* pos = (Magnum::Vector3*)data;
            return pos[i];
        };
        _crossSections[id].nematic->scale *= next;
        _crossSections[id].nematic->Init(positions.get(), access, size_reduced);
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
        if (ImGui::BeginMenu("Tools")) {

            if (ImGui::MenuItem("Isosurfaces")) {
                _widget.showPreimageSettings = true;
            }

            if (ImGui::MenuItem("Vortex Knot")) {
                _widget.showVortexKnotSettings = true;
            }

            if (ImGui::MenuItem("Nonlinear Imaging")) {
                _widget.showNonlinearSettings = true;
            }

            if (ImGui::MenuItem("Z-Profile")) {
                _widget.showZProfileWindow = true;
            }

            if (ImGui::MenuItem("Knot Interaction")) {
                _widget.knot_interaction_handle.showWindow = true;
            }

            if (ImGui::MenuItem("Multi plane")) {
                _widget.multiplane_window = true;
            }

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Settings")) {

            if (ImGui::MenuItem("POM")) {
                // Toggle POM window
                _widget.showPOMSettings = true;
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

        if (ImGui::Button("Create Hopf Density Isosurface")) {
            _hopfionDensityIsosurface = HopfDensityIsosurface{};
        }

        // Quaternion preimages
        ImGui::PushItemWidth(200.0f);
        ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "Pion field preimages");
        ImGui::SliderFloat3("pi", &_widget.pionComponents[0], 0, 1.0f);
        ImGui::PopItemWidth();

        if (ImGui::Button("Add pion preimage")) {
            _pionPreimage = PionPreimage{_widget.pionComponents, _widget.isoLevel};
        }

        ImGui::SameLine();

        if (ImGui::Button("Remove pion preimage")) {
            if (_pionPreimage)
                _pionPreimage->draw = false;
            _pionPreimage = {};
        }

        if (_pionPreimage) {
            ImGui::SameLine();
            ImGui::PushItemWidth(50.0f);
            ImGui::SliderFloat("Alpha", &_pionPreimage->alpha, 0.f, 1.f);
            ImGui::PopItemWidth();
            ImGui::SameLine();
            ImGui::Checkbox("Draw", &_pionPreimage->draw);
        }
       
        
        if (ImGui::Button("Generate Isosurfaces") || _widget.generateKnots) {
            generateIsosurface();
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

            }

            _widget.updateImage = true;
            _widget.savePOM = false;
        }
    }
}

void Sandbox::handleVortexKnotWindow() {
    if (_widget.showVortexKnotSettings) {
        Dataset* data = (Dataset*)(_solver->GetDataPtr());

        ImGui::Begin("Knot Search Tool", &_widget.showVortexKnotSettings);

        
        ImGui::PushItemWidth(120.f);

        if (ImGui::Button("Find Components")) {
            findVortexKnotComponents2();
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

        if (_widget.generateKnots) {
            generateIsosurface();
        }

        ImGui::Checkbox("Draw Surfaces", &_widget.drawSurfaces);

        ImGui::End();

    }
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
            
            if (ImPlot::BeginPlot("Free Energy v. Iterations", "Iterations", "Free energy (Kp)")) {
                ImPlot::SetNextMarkerStyle(ImPlotMarker_None);

                // Reset points
                if (_widget.updateImageFromLoad) points = MAX_GRAPH_POINTS;

                ImPlot::PlotLine("F", &_widget.series_x_axis_vec[0], &_widget.energy_series_vec[0], points, 0, sizeof(LC::precision_scalar));
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
                //std::string c1 = std::to_string(com(0)) + ", " + std::to_string(com(1)) + ", " + std::to_string(com(2));
                //std::string tstr = std::to_string(translation(0)) + ", " + std::to_string(translation(1)) + ", " + std::to_string(translation(2));
                //std::cout << c1 << std::endl;
                displacement_vectors.emplace_back(com);
            }
            else { // Invalid data point return empty vector
                return heliknotons_found;
            }
        }

    }

    // <missile>
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

        static float phase_offset = 0.f;

        ImGui::InputFloat("Phase offset", &phase_offset);
        

        ImGui::RadioButton("Toron", &_widget.hopfion_type, 2);
        ImGui::SameLine();
        ImGui::RadioButton("Twistion (T1B2)", &_widget.hopfion_type, 3);
        ImGui::SameLine();
        ImGui::RadioButton("Twistion (T2B2)", &_widget.hopfion_type, 4);

        ImGui::RadioButton("CF3", &_widget.hopfion_type, 5);
        ImGui::SameLine();
        ImGui::RadioButton("CF1", &_widget.hopfion_type, 6);
        ImGui::SameLine();
        ImGui::RadioButton("CF1 (L)", &_widget.hopfion_type, 7);

        int Q = _widget.topological_charge;
        int npp = _widget.npp;
        float v0 = _widget.voltage;

        if (ImGui::Button("Generate")) {

            std::array<int, 3> V = { _widget.celldims[0] * npp, _widget.celldims[1] * npp, _widget.celldims[2] * npp };

            (*data).Cell(_widget.celldims[0], _widget.celldims[1], _widget.celldims[2]);

            Dataset::Config cfg;

            if (_widget.hopfion_type == 0) cfg = Dataset::Heliknoton(Q, V, data->cell_dims,1.,1.,{0.f,0.f,0.f},phase_offset);
            else if (_widget.hopfion_type == 1) cfg = Dataset::Hopfion(Q);
            else if (_widget.hopfion_type == 2) cfg = Dataset::Toron(data->cell_dims, Q);
            else if (_widget.hopfion_type == 3) cfg = Dataset::Twistion_T1B2(data->cell_dims);
            else if (_widget.hopfion_type == 4) cfg = Dataset::Twistion_T2B2(data->cell_dims);
            else if (_widget.hopfion_type == 5) cfg = Dataset::CF3(data->cell_dims[0]);
            else if (_widget.hopfion_type == 6) cfg = Dataset::CF1(data->cell_dims);
            else if (_widget.hopfion_type == 7) cfg = Dataset::CF1_Loop(data->cell_dims);
            else cfg = Dataset::Toron(data->cell_dims, Q);

            data->chirality = _widget.chirality;

            (*data).Voxels(V[0], V[1], V[2])
                .Boundaries(_widget.boundaries[0], _widget.boundaries[1], _widget.boundaries[2]);
            (*data).Configuration(cfg);

            // Reset data
            (*data).numIterations = 0;

            _widget.energy_series = std::list<LC::precision_scalar>{};
            _widget.series_x_axis = std::list<LC::precision_scalar>{};
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
            if (_widget.hopfion_type == 0) cfg = Dataset::Heliknoton(Q,data->voxels,data->cell_dims);
            else if (_widget.hopfion_type == 1) cfg = Dataset::Hopfion(Q);
            else cfg = Dataset::Toron(data->cell_dims, Q);

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

            // Normalize director field
            for (std::size_t i = 0; i < numDirectors; i++) {
                LC::scalar norm = 0.;
                for (int d = 0; d < 3; d++)
                    norm += data->directors[i + d * numDirectors] * data->directors[i + d * numDirectors];
                if (norm != 0.) {
                    norm = sqrt(norm);
                    for (int d = 0; d < 3; d++)
                        data->directors[i + d * numDirectors] /= norm;
                }
                else {
                    // Director averaged to zero -> just choose n = zhat
                    data->directors[i + 2 * numDirectors] = 1.;
                }

            }

            // Set energy again
            data->en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[numDirectors]);

            solver->Normalize();
            initVisuals();
            updateColor();
        }

        if (ImGui::CollapsingHeader("Tilt and Lehman")) {

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

                            Eigen::Quaternion<LC::scalar> vpos(0., i - (data->voxels[0] - 1) / 2., j - (data->voxels[1] - 1) / 2., k - (data->voxels[2] - 1) / 2.);

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

                            auto n = helical(pos[0], pos[1], pos[2] - _widget.helicalLayerOffset);

                            for (int d = 0; d < 3; d++) {
                                nn(i, j, k, d) = n[d];
                            }
                        }
            }

            ImGui::SameLine();
            ImGui::PushItemWidth(50.f);
            ImGui::InputFloat("Layer offset (p)", &_widget.helicalLayerOffset);
            ImGui::PopItemWidth();
            ImGui::SameLine();

            static bool s_uniformflip = false;

            if (ImGui::Button("Uniform field")) {
                FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

                std::array<float, 3> pos;

                for (int i = _widget.shrink_interval_begin[0] - 1; i < _widget.shrink_interval_end[0]; i++)
                    for (int j = _widget.shrink_interval_begin[1] - 1; j < _widget.shrink_interval_end[1]; j++)
                        for (int k = _widget.shrink_interval_begin[2] - 1; k < _widget.shrink_interval_end[2]; k++) {
                            nn(i, j, k, 0) = 0;
                            nn(i, j, k, 1) = 0;
                            if (s_uniformflip)
                                nn(i, j, k, 2) = -1.0;
                            else
                                nn(i, j, k, 2) = 1.0;
                        }
            }
            ImGui::SameLine();
            ImGui::Checkbox("Invert", &s_uniformflip);



            // Resets/initializes lehman cluster
            if (ImGui::Button("Fill Lehman cluster")) {
                _lehman_xsection.fill(data->voxels, data->cell_dims, _widget.lehman_upsample);
                uint vol = data->voxels[0] * data->voxels[1] * data->voxels[2];
                for (uint i = 0; i < _lehman_xsection.size(); i++) {
                    uint id = _lehman_xsection.index(i);
                    data->directors[id] = _lehman_xsection.nodes[i].director.x();
                    data->directors[id + vol] = _lehman_xsection.nodes[i].director.y();
                    data->directors[id + 2 * vol] = _lehman_xsection.nodes[i].director.z();
                }
            }
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
                z_translation = true;
            }

            ImGui::SameLine();

            int sgn = 1;

            if (ImGui::Button("Down")) {
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


            // Prototyping
            if (ImGui::Button("Lehman cluster line")) {
                // Default field configuration: ng = (-sin(qz), cos(qz), 0) so ng(z=0) = (0,1,0)
                // Want the Lehman cluster to be in the xz plane
                auto helical = LC::Math::Planar(2, 1);

                FOFDSolver::Tensor4 nn(data->directors.get(), data->voxels[0], data->voxels[1], data->voxels[2], 3);

                LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
                LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);

                // Input:
                // r0: Position of the defect
                // y0: Y plane where the defect is placed
                // sgn: Which side the defect is placed (+ right, - left)
                auto LehmanCluster = [&](Eigen::Vector3d r0, int y0, int sgn) {

                    std::array<LC::scalar, 3> n0 = helical(r0.x(), r0.y(), r0.z());

                    for (int x = 0; x < data->voxels[0]; x++) {
                        for (int z = 0; z < data->voxels[2]; z++) {
                            // Get current position
                            Eigen::Vector3d pos(
                                -data->cell_dims[0] * 0.5 + x * dx,
                                0.,
                                -data->cell_dims[2] * 0.5 + z * dz
                            );

                            // Position relative to the point defect
                            Eigen::Vector3d rprime = pos - r0;
                            LC::scalar rprime_len = rprime.norm();

                            // If within half a pitch from the point defect and to the right of the defect
                            if (rprime_len > 0. && rprime_len <= 0.5 && sgn * rprime.x() > 0.) {
                                Eigen::Quaterniond yhat(0., n0[0], n0[1], n0[2]); // Initial director orientation at center of defect
                                // Angle of position relative to point defect in the defect plane
                                LC::scalar phi = atan2(rprime.z(), rprime.x());

                                // Radial rotation quaternion
                                LC::scalar rprime_len = rprime.norm();
                                LC::scalar theta = 2. * M_PI * rprime_len;
                                LC::scalar ct, st;
                                ct = cos(0.5 * theta);
                                st = sin(0.5 * theta);
                                Eigen::Quaterniond rot_quat;
                                rot_quat.w() = ct;
                                rot_quat.x() = st * rprime.x() / rprime_len;
                                rot_quat.y() = st * rprime.y() / rprime_len;
                                rot_quat.z() = st * rprime.z() / rprime_len;

                                // Apply the quaternion to phihat
                                Eigen::Quaterniond n = rot_quat * yhat * rot_quat.conjugate();
                                nn(x, y0, z, 0) = n.x();
                                nn(x, y0, z, 1) = n.y();
                                nn(x, y0, z, 2) = n.z();

                            }

                        }

                    }
                };

                for (int y = 0.25 * data->voxels[1]; y < 0.75 * data->voxels[1]; y++) {
                    // Lehman cluster at x = -0.5
                    LehmanCluster({ -0.5, 0., 0. }, y, 1);
                    // Lehman cluster at x = 0.5
                    LehmanCluster({ 0.5, 0., 0. }, y, -1);
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
                        for (const auto& k : loop) {

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
        }

        // 0 -> spherical interaction energy
        // 1 -> radial interaction energy
        static int radioInteractionType = 0;
        static float a_axis = 1.f;
        static float b_axis = 1.f;
        static float c_axis = 1.1f;
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


        ImGui::PushItemWidth(100.f);
        if (radioInteractionType == 0) {
            ImGui::InputFloat("Separation Distance (x)", &_widget.separationDistancex);
            ImGui::InputFloat("Separation Distance (y)", &_widget.separationDistancey);
            ImGui::InputFloat("Separation Distance (z)", &_widget.separationDistancez);
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
        ImGui::Checkbox("Force cell size", &_widget.forceInteractionCellSize);
        ImGui::SameLine();
        ImGui::Checkbox("Single Heliknoton", &_widget.singleInteractionHeliknoton);
        ImGui::SameLine();
        ImGui::PushItemWidth(150.f);
        ImGui::InputFloat3("Cell size", &_widget.interactionCellSize[0]);
        ImGui::PopItemWidth();
        ImGui::InputFloat("E-field", &EField);
        ImGui::InputInt("Relaxation iterations", &_widget.interactionIterations);
        ImGui::InputInt("Samples per point", &Navg);
        ImGui::PopItemWidth();

        // E-field

        const size_t interaction_fname_buffer_size = 128;
        static char interaction_fname[interaction_fname_buffer_size] = "interaction.bin";

        ImGui::InputText("Interaction file name", interaction_fname, interaction_fname_buffer_size);

        if (ImGui::Button("Heliknoton interaction landscape")) {

            // Get all points belonging to the heliknoton
            // ============================================
            unsigned int vol = data->voxels[0] * data->voxels[1] * data->voxels[2];
            unsigned int bulk_vol = (data->voxels[0]-2) * (data->voxels[1]-2) * (data->voxels[2]-2);
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

            Position center_of_mass(0.,0.,0);

            // Extract the heliknoton
            for (int i = 0; i < data->voxels[0]; i++) {
                for (int j = 0; j < data->voxels[1]; j++) {
                    for (int k = 0; k < data->voxels[2]; k++) {
                        // Diffmag for helical axis
                        unsigned int idx = i + data->voxels[0] * j + data->voxels[0] * data->voxels[1] * k;

#if 0 // Use handedness tensor
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
                        if (R2 > 0.2) {
                            // Add location to the heliknoton
                            Position pos;
                            pos[0] = data->cell_dims[0] * ((LC::scalar)i / (data->voxels[0] - 1) - 0.5);
                            pos[1] = data->cell_dims[1] * ((LC::scalar)j / (data->voxels[1] - 1) - 0.5);
                            pos[2] = data->cell_dims[2] * ((LC::scalar)k / (data->voxels[2] - 1) - 0.5);

                            center_of_mass += pos;

                            Director dir;
                            dir[0] = data->directors[idx];
                            dir[1] = data->directors[idx + vol];
                            dir[2] = data->directors[idx + 2 * vol];
                            heliknoton.push_back({ pos, dir });
                        }
#else   // Use dot product relative to the helical background
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
                        if (pow(pos[0]/a_axis,2) + pow(pos[1]/b_axis,2) + pow(pos[2]/c_axis,2) <= 1)
                            heliknoton.push_back({ pos, dir });
#endif

                       
                    }
                }
            }

            if (!heliknoton.empty())
                center_of_mass = center_of_mass / heliknoton.size();
            else {
                LC_ERROR("Empty heliknoton!!!");
                return;
            }


            //LC_INFO("Heliknoton center of mass R = ({0}, {1}, {2})", center_of_mass(0), center_of_mass(1), center_of_mass(2));

            // Translate the heliknoton center of mass to the origin
            for (auto& p : heliknoton) {
                p.position = p.position - center_of_mass;
            }

            // Get a bounding box on the heliknoton size
            Position lower_bound(0.,0.,0.), upper_bound(0.,0.,0.);
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
                cdims[0] = 2. * background_padding * (upper_bound[0] - lower_bound[0]) + 0.5 * _widget.separationDistancex;
                cdims[1] = 2. * background_padding * (upper_bound[1] - lower_bound[1]) + 0.5 * _widget.separationDistancey;
                cdims[2] = 2. * background_padding * (upper_bound[2] - lower_bound[2]) + 0.5 * _widget.separationDistancez;
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

                //translation_vector[0] = translation_vector[0] / data->cell_dims[0] * 3;
                //translation_vector[1] = translation_vector[1] / data->cell_dims[1] * 3;
                //translation_vector[2] = translation_vector[2] / data->cell_dims[2] * 3;

                //auto heliknoton_cfg = FOFDSolver::Dataset::Heliknoton(1, 3. / data->cell_dims[2], 1., { 0.0,0.0,0.0 }, false);
                LC::scalar omega0 = 2. * M_PI * data->chirality * translation_vector[2] + _widget.interactionOmegaOffset * M_PI / 180.f;
                auto heli_vfield = LC::Math::Heliknoton(1, data->cell_dims, 1., 1., { 0.,0.,0. }, omega0, false);
                //std::vector<uint32_t> update_list;

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
                // Standard deviation of the energy
                LC::scalar x_1f; LC::scalar y_1f; LC::scalar z_1f;
                LC::scalar x_2f; LC::scalar y_2f; LC::scalar z_2f;
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

                        LC::scalar rx = 0.5 * _widget.separationDistancex;
                        LC::scalar ry = 0.5 * _widget.separationDistancey;
                        LC::scalar rz = 0.5 * _widget.separationDistancez;

                        Eigen::Vector3d disp;
                        disp[0] = rx * sin(theta) * cos(phi);
                        disp[1] = ry * sin(theta) * sin(phi);
                        disp[2] = rz * cos(theta);

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
                    displacement[0] = rr * sin(theta0*M_PI/180.f) * cos(phi0*M_PI/180.f);
                    displacement[1] = rr * sin(theta0*M_PI/180.f) * sin(phi0*M_PI/180.f);
                    displacement[2] = rr * cos(theta0*M_PI/180.f);

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
            solver->SetVoltage(EField* data->cell_dims[2], 50);
            LC::scalar E_bg = kbT_conversion * solver->TotalEnergy();
            LC::scalar E_H;
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
                        interaction_data[count].separation *= 0.5; // Need to multiply by half to get the distance from the center
                        interaction_data[count].x_1f = 0.0; 
                        interaction_data[count].y_1f = 0.0; 
                        interaction_data[count].z_1f = 0.0; 
                        interaction_data[count].x_2f = 0.0; 
                        interaction_data[count].y_2f = 0.0; 
                        interaction_data[count++].z_2f = 0.0; 
                    }
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

                    std::vector<Eigen::Vector3d> translation_pair;
                    translation_pair.push_back(translation);
                    translation_pair.push_back(-translation);
                    
                    // owo
                    // TODO: Reimplement heliknoton code, need to compute preimage then pass to detect heliknoton preimage
                    
                    LC_INFO("Energy: {2}, Iteration {0}/{1}", count + 1, translations.size(), interaction_data[count].energy);
                    ++count;
                }
            }

            // Save data
            //std::string saveName = "D:\\dev\\lclab2\\data\\interactions\\" + std::string(interaction_fname);
            
         
            std::string saveName = "D:\\dev\\lclab2\\data\\interactions\\" + std::string(interaction_fname);
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



void Sandbox::findVortexKnotComponents(bool saveObj, const std::string &obj_name) {
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

    std::unique_ptr<float[]> chi_field, field_nn(new float[3 * vol]), valid_fieldf(new float[vol]), field_S(new float[vol]);

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

    _vortex_line.numNodes = LC::Math::ChiralityField(field_nn.get(), chi_field, data->voxels, cellf, _vortex_line.valid_field, true, _widget.knot_interaction_handle.tolerance, false);

    // Convert valid field to floats and make a list of all vortex indices
    std::vector<std::size_t> vortex_indices;
    std::set<std::size_t> vortex_indices_set;
    vortex_indices.reserve(_vortex_line.numNodes);

#define RBF_SAMPLED_VORTEX 0
#if RBF_SAMPLED_VORTEX

    std::set<std::size_t> vortex_indices_extended_set, vortex_indices_extended_subset, parsed_vortex_indices;
    std::vector<std::size_t> vortex_indices_extended, vortex_indices_extended_sublist;

    for (auto i = 0; i < vol; i++) {
        if (!_vortex_line.valid_field[i]) // Not a valid vortex line point
            vortex_indices.emplace_back(i);
    }

    int Radius = 8;
    int SubRadius = 3;
    // For each defect, increase the weight by radius Radius
    for (const auto& id : vortex_indices) {
        // Convert the index to i,j,k space
        int k = id / slice;
        int j = (id - k * slice) / data->voxels[0];
        int i = id - j * data->voxels[1] - slice * k;

        for (int x = -Radius; x <= Radius; x++)
            for (int y = -Radius; y <= Radius; y++)
                for (int z = -Radius; z <= Radius; z++) {
                    // Compute the index to modify
                    uint id_ball = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                    // Make sure defect is in bounds
                    if (id_ball < vol && id_ball >= 0) {
                        vortex_indices_extended_set.insert(id_ball);
                        if (abs(x) <= SubRadius && abs(y) <= SubRadius && abs(z) <= SubRadius)
                            vortex_indices_extended_subset.insert(id_ball);
                    }
                }
    }

    // Fill the extended index list
    vortex_indices_extended.reserve(vortex_indices_extended_set.size());
    for (auto id : vortex_indices_extended_set) {
        vortex_indices_extended.emplace_back(id);
    }
    {
        int counter = 0;
        for (const auto& id : vortex_indices_extended) {
            // Search for the extended set index in the extended subset
            auto search = vortex_indices_extended_subset.find(id);
            if (search != vortex_indices_extended_subset.end()) // Found the point in subset
                vortex_indices_extended_sublist.emplace_back(counter);
            ++counter;
        }
    }


    // Compute the position for each point in the extended knot
    std::unique_ptr<LC::scalar[]> extended_knot_pos(new LC::scalar[3 * vortex_indices_extended.size()]);

    {
        LC::scalar dx = data->cell_dims[0] / (data->voxels[0] - 1);
        LC::scalar dy = data->cell_dims[1] / (data->voxels[1] - 1);
        LC::scalar dz = data->cell_dims[2] / (data->voxels[2] - 1);
        int counter = 0;
        for (auto id : vortex_indices_extended) {
            int k = id / slice;
            int j = (id - k * slice) / data->voxels[0];
            int i = id - j * data->voxels[1] - slice * k;
            LC::scalar x = i * dx - data->cell_dims[0] * 0.5;
            LC::scalar y = j * dy - data->cell_dims[1] * 0.5;
            LC::scalar z = k * dz - data->cell_dims[2] * 0.5;

            extended_knot_pos[counter] = x;
            extended_knot_pos[counter + vortex_indices_extended.size()] = y;
            extended_knot_pos[counter + 2 * vortex_indices_extended.size()] = z;
            ++counter;
        }
    }

    // Reclassify vortex indices by interpolating and recomputing discriminant
    int knn = 60;
    LC::Math::Metric<LC::scalar> metric;
    LC::Math::StencilWeights<LC::scalar> derivative = LC::Math::StencilWeightFirstDerivatives<LC::scalar>{};
    metric.Bcs = {false, false, false};
    metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);
    std::unique_ptr<LC::Math::rbf<LC::scalar>> RBF = std::unique_ptr<LC::Math::poly_spline<LC::scalar>>(new LC::Math::poly_spline<LC::scalar>);
    std::unique_ptr<std::size_t[]> neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[vortex_indices_extended_sublist.size() * knn]);

    // Find nearest neighbors
    LC::Algorithm::knn_c(extended_knot_pos.get(), vortex_indices_extended.size(),
        &vortex_indices_extended_sublist[0], vortex_indices_extended_sublist.size(),
        metric, knn, (LC::scalar*)0, neighbors.get());

    // Compute weights
    LC_INFO("sublist/total= {0}/{1}", vortex_indices_extended_sublist.size(), vortex_indices_extended.size());

    derivative.ComputeWeights(extended_knot_pos.get(), neighbors.get(), *(RBF.get()), metric,
        vortex_indices_extended_sublist.size(), vortex_indices_extended.size(), knn);

    // Extract weights
    LC::Math::Weight<LC::scalar>* w_dx = derivative.GetWeight(LC::Math::WeightTag::x);
    LC::Math::Weight<LC::scalar>* w_dy = derivative.GetWeight(LC::Math::WeightTag::y);
    LC::Math::Weight<LC::scalar>* w_dz = derivative.GetWeight(LC::Math::WeightTag::z);

    LC::scalar nx, ny, nz;

    for (std::size_t index = 0; index < vortex_indices_extended_sublist.size(); index++) {
        // Compute derivatives
        int mapoffset = knn * index;
        LC::scalar nx100(0.), nx010(0.), nx001(0.), ny100(0.), ny010(0.), ny001(0.), nz100(0.), nz010(0.), nz001(0.);

        std::size_t nbh;
        for (int i = 0; i < knn; i++) {

            nbh = vortex_indices_extended[neighbors[i * vortex_indices_extended_sublist.size() + index]];


            nx = field_nn[nbh];
            ny = field_nn[nbh + vol];
            nz = field_nn[nbh + 2 * vol];
            //LC_INFO("dx {0}, dy {1}, dz {2}", w_dx->data[mapoffset + i], w_dy->data[mapoffset + i], w_dz->data[mapoffset + i]);

            nx100 += w_dx->data[mapoffset + i] * nx;
            ny100 += w_dx->data[mapoffset + i] * ny;
            nz100 += w_dx->data[mapoffset + i] * nz;

            nx010 += w_dy->data[mapoffset + i] * nx;
            ny010 += w_dy->data[mapoffset + i] * ny;
            nz010 += w_dy->data[mapoffset + i] * nz;

            nx001 += w_dz->data[mapoffset + i] * nx;
            ny001 += w_dz->data[mapoffset + i] * ny;
            nz001 += w_dz->data[mapoffset + i] * nz;
        }

        // Recompute the discriminant
        LC::scalar new_discr = -pow((-nx010 + ny100) * nz + nx * (-ny001 + nz010) + ny * (nx001 - nz100), 2) +
            2 * (pow(nx001 * ny - nx * ny001, 2) + 2 * (-(nx010 * ny) + nx * ny010) * (nx001 * nz - nx * nz001) + 2 * (nx100 * ny - nx * ny100) * (ny001 * nz - ny * nz001) + pow(nx010 * nz - nx * nz010, 2) +
                2 * (-(ny010 * nz) + ny * nz010) * (nx100 * nz - nx * nz100) + pow(ny100 * nz - ny * nz100, 2));

        LC_INFO("Discriminant = {0}", new_discr);

        // New point is belongs to vortex list
        if (new_discr <= _widget.knot_interaction_handle.isoValue)
            parsed_vortex_indices.insert(vortex_indices_extended[vortex_indices_extended_sublist[index]]);
    }

    LC_INFO("Number of resolved vortex line points = {0}", parsed_vortex_indices.size());

    for (auto i = 0; i < vol; i++) {
        auto search = parsed_vortex_indices.find(i);
        if (search == parsed_vortex_indices.end()) // Not a valid vortex line point
            valid_fieldf[i] = 1.f;
        else { // Vortex line
            valid_fieldf[i] = -1.f;
        }
    }
#else
    // Convert valid field to floats and make a list of all defect indices
    for (auto ii = 0; ii < vol; ii++) {
        valid_fieldf[ii] = 0.f;
        if (_vortex_line.valid_field[ii]) {

            int k = ii / slice;
            int j = (ii - k * slice) / data->voxels[0];
            int i = ii - j * data->voxels[1] - slice * k;
            LC::scalar avgcos2_chi = 0.;
            LC::scalar avgcos2_tau = 0.;
            Eigen::Vector3d chi0(chi_field[ii], chi_field[ii + vol], chi_field[ii + 2 * vol]);
            Eigen::Vector3d nn0(field_nn[ii], field_nn[ii + vol], field_nn[ii + 2 * vol]);
            Eigen::Vector3d tau0 = chi0.cross(nn0);
            // Compute the scalar order parameter with the chi field


            int points_in_bds = 0;

            for (int x = -1; x <= 1; x++) {
                for (int y = -1; y <= 1; y++) {
                    for (int z = -1; z <= 1; z++) {
                        uint ip = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                        if (ip >= 0 && ip < vol && ip != ii) {
                            points_in_bds++;
                            Eigen::Vector3d chi_xyz(chi_field[ip], chi_field[ip + vol], chi_field[ip + 2 * vol]);
                            Eigen::Vector3d nn_xyz(field_nn[ip], field_nn[ip + vol], field_nn[ip + 2 * vol]);
                            Eigen::Vector3d tau_xyz = chi_xyz.cross(nn_xyz);
                            avgcos2_chi += pow(chi0.dot(chi_xyz), 2);
                            avgcos2_tau += pow(tau0.dot(tau_xyz), 2);
                        }
                    }
                }
            }

            LC::scalar Sop_chi,Sop_lambda,Sop_tau;
            if (points_in_bds > 0) {
                Sop_chi = 0.5 * (3. * avgcos2_chi / points_in_bds - 1.);
                Sop_tau = 0.5 * (3. * avgcos2_tau / points_in_bds - 1.);
                LC::scalar S = 0.5 * (Sop_chi + Sop_tau);
                field_S[ii] = S;
                if (S < _widget.knot_interaction_handle.isoValue) {
                    //vortex_indices.emplace_back(ii);
                    //vortex_indices_set.insert(ii);
                    _vortex_line.valid_field[ii] = 0;
                }
            }

        }
        else {
            //vortex_indices.emplace_back(ii);
            field_S[ii] = 0.f;
        }
    }

// Try secondary vortex knot classification based on |grad(S)|=0 and S <= S.isoValue

    for (auto ii = 0; ii < vol; ii++) {
        int k = ii / slice;
        int j = (ii - k * slice) / data->voxels[0];
        int i = ii - j * data->voxels[1] - slice * k;

        int line = data->voxels[0];

        auto get_id = [slice, line](int x, int y, int z) {
            return x + line * y + slice * z;
        };

        // Compute S derivatives
        float Sd2 = 0.f;
        if (i > 1 && j > 1 && k > 1 && i < data->voxels[0] - 2 && j < data->voxels[1] - 2 && k < data->voxels[2] - 2) {
            float S200 = (field_S[get_id(i + 1, j, k)] + field_S[get_id(i - 1, j, k)] - 2. * field_S[ii]);
            float S020 = (field_S[get_id(i, j + 1, k)] + field_S[get_id(i, j - 1, k)] - 2. * field_S[ii]);
            float S002 = (field_S[get_id(i, j, k + 1)] + field_S[get_id(i, j, k - 1)] - 2. * field_S[ii]);
            Sd2 = S200 + S020 + S002;
        }


        // Classify the current index ii
        if (Sd2 > 0.5 && field_S[ii] < _widget.knot_interaction_handle.isoValue) {
            vortex_indices.emplace_back(ii);
            vortex_indices_set.insert(ii);
        }

    }

    int Radius = _widget.knot_interaction_handle.saturation;
    int RadiusSq = Radius * Radius;
    // For each defect point, enhance contrast by expanding outwards

    for (const auto& id : vortex_indices) {
        // Convert the index to i,j,k space
        int k = id / slice;
        int j = (id - k * slice) / data->voxels[0];
        int i = id - j * data->voxels[0] - slice * k;

        for (int x = -Radius; x <= Radius; x++)
            for (int y = -Radius; y <= Radius; y++)
                for (int z = -Radius; z <= Radius; z++) {
                    // Compute the index to modify
                    uint id_ball = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                    // Make sure defect is in bounds
                    if (id_ball < vol && id_ball >= 0) {
                        //_vortex_line.valid_field[id_ball] = 0;
                        float r2 = x * x + y * y + z * z;

                        if (id_ball == id)
                            valid_fieldf[id_ball] += -1.f;
                        else
                            valid_fieldf[id_ball] += -1.f / r2;
                    }
                }

    }

    // Recount newly classified nodes
    _vortex_line.numNodes = 0;

    for (auto ii = 0; ii < vol; ii++) {
        valid_fieldf[ii] += Radius;
        if (valid_fieldf[ii] < 0.f) {
            _vortex_line.valid_field[ii] = 0;
            _vortex_line.numNodes++;
        }
    }

#endif

    std::map<uint, Face> unvisited_list, unvisited_list_local;
    // Key: vertex, value: Triangle
    std::multimap<uint, uint> mmap;

    // Make the initial color white to distinguish components that were too small
    std::array<float, 4> col = { 1.f, 1.f, 1.f, 1.f };

    //LC::Math::Isosurface<short*, short> gen2;
    //gen2.GenerateSurface(_vortex_line.valid_field.get(), 0, data->voxels, dr, col);

    // Use extended marching cubes to generate the isosurface
    LC::ExtendedMC::MarchingCubes mc;

    mc.clean_all();
    mc.set_resolution(data->voxels[0], data->voxels[1], data->voxels[2]);
    mc.init_all();

    mc.set_ext_data(valid_fieldf.get());
    mc.run();

    // Rescale positions
    {
        float dx = data->cell_dims[0] / (data->voxels[0] - 1);
        float dy = data->cell_dims[1] / (data->voxels[1] - 1);
        float dz = data->cell_dims[2] / (data->voxels[2] - 1);
        float xmin = -0.5 * data->cell_dims[0];
        float ymin = -0.5 * data->cell_dims[1];
        float zmin = -0.5 * data->cell_dims[2];

        for (uint i = 0; i < mc.nverts(); ++i) {
            LC::ExtendedMC::Vertex& v = mc.vertices()[i];
            v.x = dx * v.x + xmin;
            v.y = dy * v.y + ymin;
            v.z = dz * v.z + zmin;

            // Normalize normals
            float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
            if (nrm != 0)
            {
                nrm = 1.0 / sqrt(nrm);
                v.nx *= nrm;
                v.ny *= nrm;
                v.nz *= nrm;
            }
        }
    }

    if (mc.nverts() && mc.ntrigs()) {
        unsigned int nVert = mc.nverts();
        unsigned int nInd = mc.ntrigs() * 3;
        LC::Math::IsoVertex* verts = new LC::Math::IsoVertex[nVert];
        unsigned int* indices = new unsigned int[nInd];

        LC::ExtendedMC::Vertex* temp_verts = mc.vertices();
        int* ind_temp = (int*)mc.triangles();

        // Feed data to indices,verts

        for (int i = 0; i < nInd; i++)
            indices[i] = ind_temp[i];

        for (int i = 0; i < nVert; i++) {
            verts[i].position[0] = temp_verts[i].x;
            verts[i].position[1] = temp_verts[i].y;
            verts[i].position[2] = temp_verts[i].z;
            verts[i].normal[0] = temp_verts[i].nx;
            verts[i].normal[1] = temp_verts[i].ny;
            verts[i].normal[2] = temp_verts[i].nz;

        }

        mc.reset_mesh();

        std::vector<MeshLib::PNCVertex<float>> vertices;
        std::vector<MeshLib::Triangle> triangles;
        {

            std::array<float, 4> randcol = { 1.f,1.f,1.f,1.f };
            auto smoothing_result = TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(), data->voxels, randcol, false, false, 0, saveObj, obj_name);

            vertices = std::get<0>(smoothing_result);
            triangles = std::get<1>(smoothing_result);

            // Pass vertices and triangles back to verts and indices
            delete[] verts;
            delete[] indices;

            nVert = vertices.size();
            nInd = triangles.size() * 3;
            
            verts = (LC::Math::IsoVertex*)&vertices[0];
            indices = (unsigned int*)&triangles[0][0];
        }


        unsigned int nTriangles = nInd / 3;

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
        try {
            _vortex_line.components = find_all_components_graph(unvisited_list_local, nVert, _widget.minComponentSize);
        }
        catch (const std::runtime_error& e) {
            LC_INFO("Error: {0}", e.what());
        }

        _vortex_line.vertices = vertices;
        // Color vertices by their component
        int col_i = 0;
        float Ddeg = 360.f / (float)_vortex_line.components.size();

        // Split apart vertex data
        for (auto& component : _vortex_line.components) {


            Color4 color = Color3::fromHsv({ Deg(Ddeg * col_i++), 1.f, 1.f });
            _vortexKnot.push_back(VortexKnot{});
            VortexKnot* knot;
            {
                // Go to the newly added vortex knot
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
                // possibilities
                // 1. Auto color by component
                // 2. Director color
                bool choice = 1;
                if (choice) {
                    // Use the vertex position to determine the director location
                    int x = (v.position[0] / cellf[0] + 0.5) * (data->voxels[0] - 1);
                    int y = (v.position[1] / cellf[1] + 0.5) * (data->voxels[1] - 1);
                    int z = (v.position[2] / cellf[2] + 0.5) * (data->voxels[2] - 1);
                    uint id = x + data->voxels[0] * y + slice * z;
                    float theta = acos(data->directors[id + 2 * vol]);
                    float phi = atan2(data->directors[id + vol], data->directors[id]);
                    auto dcol = LC::Imaging::Colors::RungeSphere(theta, phi);
                    v.color[0] = dcol[0];
                    v.color[1] = dcol[1];
                    v.color[2] = dcol[2];
                    v.color[3] = color[4];
                }
                else {
                    for (int k = 0; k < 4; k++)
                        v.color[k] = color[k];
                }

                vertices.insert({ ci, v });
            }

            std::map<uint, uint> unique_triangles, vertmap;

            // Parse vertices
            std::vector<LC::Surface::Vertex> vertlist;
            vertlist.reserve(vertices.size());
            
            uint ctr = 0;
            for (auto it = vertices.begin(); it != vertices.end(); it++) {
                vertmap.insert({ it->first, ctr++ });
                vertlist.emplace_back(it->second);
            }

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
            uint count = 0;
            for (const auto& face : unique_triangles) {
                uint tri = face.first;
                
                for (int d = 0; d < 3; d++) {
                    // Vertex index
                    uint ii = indices[3 * tri + d];

                    auto viter = vertmap.find(ii);
                    // found
                    if (viter != vertmap.end()) {
                        indlist.push_back(viter->second);
                    }

                }
            }

            // Initialize the mesh
            knot->surface.Init(&vertlist[0], vertlist.size(), &indlist[0], indlist.size(), _widget.preimage_translate);

            // Create the mesh
            knot->mesh = knot->surface.Mesh();

            // Add default vortex line to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_vortexManipulator, _phongShader, *(knot->mesh), knot->draw, knot->cullFaces, _transparentNormalDrawables };

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
        //delete[] verts;
        //delete[] indices;

        LC_INFO("Discovered {0} component(s)", _vortex_line.components.size());
        for (auto& c : _vortex_line.components)
            LC_INFO("- Component size = {0}", c.size());

    }

}

void Sandbox::findVortexKnotComponents2(bool saveObj, const std::string& obj_name) {
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

    std::unique_ptr<float[]> chi_field,
        field_nn(new float[3 * vol]),
        valid_fieldf(new float[vol]),
        field_Q(new float[5 * vol]),
        field_S(new float[vol]);

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

    _vortex_line.numNodes = LC::Math::ChiralityField(field_nn.get(), chi_field, data->voxels, cellf, _vortex_line.valid_field, true, _widget.knot_interaction_handle.tolerance, false);

    // Convert valid field to floats and make a list of all vortex indices
    std::vector<std::size_t> vortex_indices;
    std::set<std::size_t> vortex_indices_set;
    vortex_indices.reserve(_vortex_line.numNodes);

    // Convert valid field to floats and make a list of all defect indices
    for (auto ii = 0; ii < vol; ii++) {
        valid_fieldf[ii] = 0.f;
        if (_vortex_line.valid_field[ii]) {

            int k = ii / slice;
            int j = (ii - k * slice) / data->voxels[0];
            int i = ii - j * data->voxels[1] - slice * k;
            LC::scalar avgcos2_chi = 0.;
            LC::scalar avgcos2_tau = 0.;
            Eigen::Vector3d chi0(chi_field[ii], chi_field[ii + vol], chi_field[ii + 2 * vol]);
            Eigen::Vector3d nn0(field_nn[ii], field_nn[ii + vol], field_nn[ii + 2 * vol]);
            Eigen::Vector3d tau0 = chi0.cross(nn0);
            // Compute the scalar order parameter with the chi field


            int points_in_bds = 0;

            for (int x = -1; x <= 1; x++) {
                for (int y = -1; y <= 1; y++) {
                    for (int z = -1; z <= 1; z++) {
                        uint ip = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                        if (ip >= 0 && ip < vol && ip != ii) {
                            points_in_bds++;
                            Eigen::Vector3d chi_xyz(chi_field[ip], chi_field[ip + vol], chi_field[ip + 2 * vol]);
                            Eigen::Vector3d nn_xyz(field_nn[ip], field_nn[ip + vol], field_nn[ip + 2 * vol]);
                            Eigen::Vector3d tau_xyz = chi_xyz.cross(nn_xyz);
                            avgcos2_chi += pow(chi0.dot(chi_xyz), 2);
                            avgcos2_tau += pow(tau0.dot(tau_xyz), 2);
                        }
                    }
                }
            }

            LC::scalar Sop_chi, Sop_lambda, Sop_tau;
            if (points_in_bds > 0) {
                Sop_chi = 0.5 * (3. * avgcos2_chi / points_in_bds - 1.);
                Sop_tau = 0.5 * (3. * avgcos2_tau / points_in_bds - 1.);
                LC::scalar S = 0.5 * (Sop_chi + Sop_tau);
                field_S[ii] = S;
                if (S < _widget.knot_interaction_handle.isoValue) {
                    _vortex_line.valid_field[ii] = 0;
                }
            }

        }
        else {
            //vortex_indices.emplace_back(ii);
            field_S[ii] = 0.f;
        }
    }

    // Try secondary vortex knot classification based on |grad(S)|=0 and S <= S.isoValue

    for (auto ii = 0; ii < vol; ii++) {
        int k = ii / slice;
        int j = (ii - k * slice) / data->voxels[0];
        int i = ii - j * data->voxels[1] - slice * k;

        int line = data->voxels[0];

        auto get_id = [slice, line](int x, int y, int z) {
            return x + line * y + slice * z;
        };

        // Compute S derivatives
        float Sd2 = 0.f;
        if (i > 1 && j > 1 && k > 1 && i < data->voxels[0] - 2 && j < data->voxels[1] - 2 && k < data->voxels[2] - 2) {
            float S200 = (field_S[get_id(i + 1, j, k)] + field_S[get_id(i - 1, j, k)] - 2. * field_S[ii]);
            float S020 = (field_S[get_id(i, j + 1, k)] + field_S[get_id(i, j - 1, k)] - 2. * field_S[ii]);
            float S002 = (field_S[get_id(i, j, k + 1)] + field_S[get_id(i, j, k - 1)] - 2. * field_S[ii]);
            Sd2 = S200 + S020 + S002;
        }


        // Classify the current index ii
        if (field_S[ii] < _widget.knot_interaction_handle.isoValue) {
            vortex_indices.emplace_back(ii);
            vortex_indices_set.insert(ii);
        }

    }

    // Classify only the raw defect core
    for (const auto& id : vortex_indices) {
        valid_fieldf[id] = -1.f;
    }

    // Recount newly classified nodes
    _vortex_line.numNodes = 0;

    for (auto ii = 0; ii < vol; ii++) {
        if (valid_fieldf[ii] < 0.f) {
            _vortex_line.valid_field[ii] = 0;
            _vortex_line.numNodes++;
        }
    }

    std::map<uint, Face> unvisited_list_local;
    // Key: vertex, value: Triangle
    std::multimap<uint, uint> mmap;
    // Make the initial color white to distinguish components that were too small
    std::array<float, 4> col_white = { 1.f, 1.f, 1.f, 1.f };

    auto extended_to_iso = [](const LC::ExtendedMC::Vertex& v) {
        LC::Math::IsoVertex viso;
        viso.position[0] = v.x; viso.position[1] = v.y; viso.position[2] = v.z;
        viso.normal[0] = v.nx; viso.normal[1] = v.ny; viso.normal[2] = v.nz;

        return viso;
    };

    // Use extended marching cubes to generate the isosurface
    LC::ExtendedMC::MarchingCubes mc, mc2;

    mc.clean_all();
    mc.set_resolution(data->voxels[0], data->voxels[1], data->voxels[2]);
    mc.init_all();

    mc.set_ext_data(valid_fieldf.get());
    mc.run();

    // Smooth the vortex knot data
    std::vector<MeshLib::PNCVertex<float>> smooth0_vertices;
    std::vector<MeshLib::Triangle> smooth0_trigs;
    {
        // Re-interleave data to IsoVertex format
        std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc.nverts()]);
        LC::ExtendedMC::Vertex* iso_v_ptr = mc.vertices();
        uint* iso_ind_ptr = (uint*)mc.triangles();

        for (int vi = 0; vi < mc.nverts(); vi++) {
            iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
        }

        LC::ExtendedMC::Triangle* tri_data = mc.triangles();

        // Random color
        std::array<float, 4> col_arr;
        for (int d = 0; d < 4; d++)
            col_arr[d] = col_white[d];

        // Apply Taubin smoothing to the pre-sampled vertices
        float lambda = 0.5f;
        float mu = -0.34f;
        int smoothIterations = 80;

        auto smoothing_result = TaubinSmoothingGeneral(lambda,mu,smoothIterations,iso_verts.get(), mc.nverts(), iso_ind_ptr,
            mc.ntrigs() * 3, 0, data->voxels, col_arr, 0, 0, 0, 0);

        smooth0_vertices = std::get<0>(smoothing_result);
        smooth0_trigs = std::get<1>(smoothing_result);

        // Clean the mc data
        mc.reset_mesh();
    }

    // Rescale positions
    {
        float dx = data->cell_dims[0] / (data->voxels[0] - 1);
        float dy = data->cell_dims[1] / (data->voxels[1] - 1);
        float dz = data->cell_dims[2] / (data->voxels[2] - 1);
        float xmin = -0.5 * data->cell_dims[0];
        float ymin = -0.5 * data->cell_dims[1];
        float zmin = -0.5 * data->cell_dims[2];

        for (uint i = 0; i < smooth0_vertices.size(); ++i) {
            MeshLib::PNCVertex<float>& v = smooth0_vertices[i];
            
            v.position[0] = dx * v.position[0] + xmin;
            v.position[1] = dy * v.position[1] + ymin;
            v.position[2] = dz * v.position[2] + zmin;

            // Normalize normals
            float nrm = v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] + v.normal[2] * v.normal[2];
            if (nrm != 0)
            {
                nrm = 1.0 / sqrt(nrm);
                v.normal[0] *= nrm;
                v.normal[1] *= nrm;
                v.normal[2] *= nrm;
            }
        }
    }

    std::array<int, 3> vox_interp;
    for (int d = 0; d < 3; d++)
        vox_interp[d] = data->voxels[d] * _widget.knot_interaction_handle.upsampling_multiplier; // Can play around with this parameter
    std::array<float, 3> dr_interp;
    for (int d = 0; d < 3; d++) {
        dr_interp[d] = cellf[d] / (vox_interp[d] - 1);
    }
    uint vol_interp = vox_interp[0] * vox_interp[1] * vox_interp[2];
    std::unique_ptr<float[]> sample_grid(new float[vol_interp]), field_Ssb(new float[vol_interp]);
    std::set<uint> vortex_indices_upsampled;

    if (smooth0_vertices.size() && smooth0_trigs.size()) {
        unsigned int nVert = smooth0_vertices.size();
        unsigned int nInd = smooth0_trigs.size() * 3;
        unsigned int nTriangles = nInd / 3;

        // Feed data to indices,verts

        //mc.reset_mesh();

        // Note that verts currently refers to verts and not global indices
        // Therefore, convert indices data to global index space through vertex positions
        for (int tri = 0; tri < nTriangles; tri++) {
            uint i1 = smooth0_trigs[tri][0];
            uint i2 = smooth0_trigs[tri][1];
            uint i3 = smooth0_trigs[tri][2];

            mmap.insert({ i1, tri });
            mmap.insert({ i2, tri });
            mmap.insert({ i3, tri });

            unvisited_list_local.insert({ tri, Face(i1, i2, i3) });
        }

        // Find components through the vertices
        try {
            _vortex_line.components = find_all_components_graph(unvisited_list_local, nVert, 15);
        }
        catch (const std::runtime_error& e) {
            LC_INFO("Error: {0}", e.what());
        }

        //_vortex_line.vertices = vertices;
        // Color vertices by their component

        int col_i = 0;
        float Ddeg = 360.f / (float)_vortex_line.components.size();

        mc2.clean_all();
        mc2.set_resolution(vox_interp[0], vox_interp[1], vox_interp[2]);
        mc2.set_ext_data(sample_grid.get());
        mc2.init_all();

        int Radius = _widget.knot_interaction_handle.saturation;


        // Split apart vertex data
        
        for (auto& component : _vortex_line.components) {

            // Initialize the sample grid
            for (int ii = 0; ii < vol_interp; ii++)
                sample_grid[ii] = 0.1f;

            // Sample the component using the indices
            for (const auto& ci : component) {
                auto v_ci = smooth0_vertices[ci];
                
                // Convert the vertex into an index on the new grid
                uint i = (v_ci.position[0] + 0.5f * cellf[0]) / cellf[0] * (vox_interp[0] - 1);
                uint j = (v_ci.position[1] + 0.5f * cellf[1]) / cellf[1] * (vox_interp[1] - 1);
                uint k = (v_ci.position[2] + 0.5f * cellf[2]) / cellf[2] * (vox_interp[2] - 1);

                // Calculate the new index
                uint id = i + j * vox_interp[0] + k * vox_interp[0] * vox_interp[1];

                // Apply some sampling technique to the vortex knot at this position
                for (int x = -Radius; x <= Radius; x++)
                    for (int y = -Radius; y <= Radius; y++)
                        for (int z = -Radius; z <= Radius; z++) {
                            // Compute the index to modify
                            uint id_ball = (i + x) + (j + y) * vox_interp[0] + (k + z) * vox_interp[0] * vox_interp[1];
                            // Make sure defect is in bounds
                            if (id_ball < vol_interp && id_ball >= 0) {
                                //_vortex_line.valid_field[id_ball] = 0;
                                float r2 = x * x + y * y + z * z;
                                if (id_ball == id)
                                    sample_grid[id_ball] += -1.f;
                                else
                                    sample_grid[id_ball] += -1.f / r2;
                            }
                        }

            }

            // Compute the new knot isosurface from sample_grid
            mc2.reset_mesh();
            mc2.run();

            // Rescale the new isosurface vertices
            {
                float xmin = -0.5 * cellf[0];
                float ymin = -0.5 * cellf[1];
                float zmin = -0.5 * cellf[2];

                for (uint i = 0; i < mc2.nverts(); ++i) {
                    LC::ExtendedMC::Vertex& v = mc2.vertices()[i];
                    v.x = dr_interp[0] * v.x + xmin;
                    v.y = dr_interp[1] * v.y + ymin;
                    v.z = dr_interp[2] * v.z + zmin;

                    // Normalize normals
                    float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
                    if (nrm != 0)
                    {
                        nrm = 1.0 / sqrt(nrm);
                        v.nx *= nrm;
                        v.ny *= nrm;
                        v.nz *= nrm;
                    }
                }
            }

            // Re-interleave data to IsoVertex format
            std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc2.nverts()]);
            LC::ExtendedMC::Vertex* iso_v_ptr = mc2.vertices();
            uint* iso_ind_ptr = (uint*)mc2.triangles();

            LC_INFO("Interp verts {0}", mc2.nverts());

            for (int vi = 0; vi < mc2.nverts(); vi++) {
                iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
            }

            LC::ExtendedMC::Triangle* tri_data = mc2.triangles();

            // Use color wheel colors for fun
            Color4 color = Color3::fromHsv({ Deg(Ddeg * col_i++), 1.f, 1.f });
            std::array<float, 4> col_arr;
            for (int d = 0; d < 4; d++)
                col_arr[d] = color[d];

            // Apply Taubin smoothing to the vertices
            float lambda = 0.5f;
            float mu = -0.34f;
            int smoothIterations = 120;
            auto smoothing_result = TaubinSmoothingGeneral(lambda,mu,smoothIterations,iso_verts.get(), mc2.nverts(), iso_ind_ptr,
                mc2.ntrigs() * 3, 0, data->voxels, col_arr, 0, 0, 0, 0);
            
            auto smooth_vertices = std::get<0>(smoothing_result);
            auto smooth_triangles = std::get<1>(smoothing_result);

            
            _vortexKnot.push_back(VortexKnot{});
            VortexKnot* knot;
            {
                // Go to the newly added vortex knot
                std::list<VortexKnot>::iterator iter;
                for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++)
                    iter = it;
                knot = &(*iter);
            }


            knot->knotColor = color;
            // Use director cols instead
            for (auto& v : smooth_vertices) {
                int x = (v.position[0] / cellf[0] + 0.5) * (data->voxels[0] - 1);
                int y = (v.position[1] / cellf[1] + 0.5) * (data->voxels[1] - 1);
                int z = (v.position[2] / cellf[2] + 0.5) * (data->voxels[2] - 1);
                int xi = (v.position[0] / cellf[0] + 0.5) * (vox_interp[0] - 1);
                int yi = (v.position[1] / cellf[1] + 0.5) * (vox_interp[1] - 1);
                int zi = (v.position[2] / cellf[2] + 0.5) * (vox_interp[2] - 1);
                uint id = x + data->voxels[0] * y + slice * z;
                uint idi = xi + vox_interp[0] * yi + vox_interp[0] * vox_interp[1] * zi;

                if (idi < vol_interp)
                    vortex_indices_upsampled.insert(idi);

                float theta = acos(data->directors[id + 2 * vol]);
                float phi = atan2(data->directors[id + vol], data->directors[id]);

                Magnum::Color3 dcol;

                dcol = LC::Imaging::Colors::RungeSphere(theta, phi);
                v.color[0] = dcol[0];
                v.color[1] = dcol[1];
                v.color[2] = dcol[2];
                v.color[3] = 1.f;
            }
            

            // Initialize the mesh
            knot->surface.Init((LC::Surface::Vertex*)&smooth_vertices[0], smooth_vertices.size(), &smooth_triangles[0][0],
                smooth_triangles.size()*3, _widget.preimage_translate);

            // Create the mesh
            knot->mesh = knot->surface.Mesh();

            // Add default vortex line to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_vortexManipulator, _phongShader, *(knot->mesh), knot->draw, knot->cullFaces, _transparentNormalDrawables };
        }

        LC_INFO("Discovered {0} component(s)", _vortex_line.components.size());
        for (auto& c : _vortex_line.components)
            LC_INFO("- Component size = {0}", c.size());

    }

    // yeet
    int Nx = data->voxels[0];
    std::set<uint> fat_vortex_list;
    auto sub2ind = [slice, Nx, vol](int i, int j, int k, int d) {
        return i + Nx * j + slice * k + vol * d;
    };
    // Compute Qij for the chi field
    for (int i = 0; i < data->voxels[0]; i++) {
        for (int j = 0; j < data->voxels[1]; j++) {
            for (int k = 0; k < data->voxels[2]; k++) {

                uint ii = i + data->voxels[0] * j + slice * k;
                // Compute the order parameter S

                Eigen::Vector3d chi0(chi_field[ii], chi_field[ii + vol], chi_field[ii + 2 * vol]);
                // Compute the scalar order parameter with the chi field

                // Use the scalar order parameter value computed above
                LC::scalar S = 1.;// field_S[ii];


                field_Q[ii] = 0.5 * S * (3. * chi0[0] * chi0[0] - 1.); // qa
                field_Q[ii + vol] = 0.5 * S * (3. * chi0[0] * chi0[1]); // qb
                field_Q[ii + 2 * vol] = 0.5 * S * (3. * chi0[0] * chi0[2]); // qc
                field_Q[ii + 3 * vol] = 0.5 * S * (3. * chi0[1] * chi0[1] - 1.); // qd
                field_Q[ii + 4 * vol] = 0.5 * S * (3. * chi0[1] * chi0[2]); // qe
            }
        }
    }

    // Fill fattened vortex excluding the knot knot itself
    int sb_rad = _widget.knot_interaction_handle.sb_saturation;
    uint slice_interp = vox_interp[0] * vox_interp[1];

    for (const auto& id : vortex_indices_upsampled) {
        // Convert the index to i,j,k space
        int k = id / slice_interp;
        int j = (id - k * slice_interp) / vox_interp[0];
        int i = id - j * vox_interp[0] - slice_interp * k;

        for (int x = -sb_rad; x <= sb_rad; x++)
            for (int y = -sb_rad; y <= sb_rad; y++)
                for (int z = -sb_rad; z <= sb_rad; z++) {
                    // Compute the index to modify
                    uint id_ball = (i + x) + (j + y) * vox_interp[0] + (k + z) * slice_interp;
                    // Make sure defect is in bounds
                    if (id_ball < vol_interp && id_ball >= 0) {
                        //_vortex_line.valid_field[id_ball] = 0;
                        // Exclude the vortex knot itself
                        if (vortex_indices_upsampled.find(id_ball) == vortex_indices_upsampled.end())
                            fat_vortex_list.insert(id_ball);
                    }
                }

    }

    // Statistical parameters
    float max_sb = 0.f;
    float min_sb = 0.f;
    float mean_sb_pos = 0.f;
    float mean_sb_neg = 0.f;
    float stdev = 0.f;
    int Npos = 0;
    int Nneg = 0;

    // Next compute the splay bend around the knot
    auto splay_bend = [&](int i, int j, int k) {
        float qa002 = (field_Q[sub2ind(i, j, k + 1, 0)] + field_Q[sub2ind(i, j, k - 1, 0)] - 2.f * field_Q[sub2ind(i, j, k, 0)]) / (dr[2] * dr[2]);
        float qa200 = (field_Q[sub2ind(i + 1, j, k, 0)] + field_Q[sub2ind(i - 1, j, k, 0)] - 2.f * field_Q[sub2ind(i, j, k, 0)]) / (dr[0] * dr[0]);
        float qd002 = (field_Q[sub2ind(i, j, k + 1, 3)] + field_Q[sub2ind(i, j, k - 1, 3)] - 2.f * field_Q[sub2ind(i, j, k, 3)]) / (dr[2] * dr[2]);
        float qd020 = (field_Q[sub2ind(i, j + 1, k, 3)] + field_Q[sub2ind(i, j - 1, k, 3)] - 2.f * field_Q[sub2ind(i, j, k, 3)]) / (dr[1] * dr[1]);
        float qb110 = (field_Q[sub2ind(i + 1, j + 1, k, 1)] + field_Q[sub2ind(i - 1, j - 1, k, 1)]
            - field_Q[sub2ind(i - 1, j + 1, k, 1)] - field_Q[sub2ind(i + 1, j - 1, k, 1)]) / (4.f * dr[0] * dr[1]);
        float qc101 = (field_Q[sub2ind(i + 1, j, k + 1, 2)] + field_Q[sub2ind(i - 1, j, k - 1, 2)]
            - field_Q[sub2ind(i - 1, j, k + 1, 2)] - field_Q[sub2ind(i + 1, j, k - 1, 2)]) / (4.f * dr[0] * dr[2]);
        float qe011 = (field_Q[sub2ind(i, j + 1, k + 1, 4)] + field_Q[sub2ind(i, j - 1, k - 1, 4)]
            - field_Q[sub2ind(i, j - 1, k + 1, 4)] - field_Q[sub2ind(i, j + 1, k - 1, 4)]) / (4.f * dr[1] * dr[2]);

        return -qa002 + qa200 + 2.f * (qb110 + qc101) - qd002 + qd020 + 2.f * qe011;
    };

    Float cijk[2][2][2];
    Float cjk[2][2];
    Float ck[2];

    for (int i = 0; i < vox_interp[0]; i++) {
        for (int j = 0; j < vox_interp[1]; j++) {
            for (int k = 0; k < vox_interp[2]; k++) {
                uint id = i + vox_interp[0] * j + slice_interp * k;
                if (i > 1 && i < vox_interp[0] - 2 && j > 1 && j < vox_interp[1] - 2 && k > 1 && k < vox_interp[2] - 2) {


                    // Found area to compute the splay bend
                    // uwu
                    if (fat_vortex_list.find(id) != fat_vortex_list.end()) {
                        // Use trilinear interpolation
                        float x = i * dr_interp[0];
                        float y = j * dr_interp[1];
                        float z = k * dr_interp[2];
                        // Bottom corner points
                        int i_0 = x / dr[0];
                        int j_0 = y / dr[1];
                        int k_0 = z / dr[2];
                        float x_0 = i_0 * dr[0];
                        float y_0 = j_0 * dr[1];
                        float z_0 = k_0 * dr[2];
                        // Interpolation coordinates
                        float x_d = (x - x_0) / dr[0];
                        float y_d = (y - y_0) / dr[1];
                        float z_d = (z - z_0) / dr[2];
                        // Fill the interpolant
                        for (int a = 0; a < 2; a++)
                            for (int b = 0; b < 2; b++)
                                for (int c = 0; c < 2; c++)
                                    cijk[a][b][c] = splay_bend(i_0 + a, j_0 + b, k_0 + c);

                        // Begin interpolation
                        for (int b = 0; b < 2; b++)
                            for (int c = 0; c < 2; c++)
                                cjk[b][c] = (1.f - x_d) * cijk[0][b][c] + x_d * cijk[1][b][c];

                        // Begin interpolation
                        for (int c = 0; c < 2; c++)
                            ck[c] = (1.f - y_d) * cjk[0][c] + y_d * cjk[1][c];

                        field_Ssb[id] = (1.f - z_d) * ck[0] + z_d * ck[1];
                    }
                    else // splay bend not close enough to the vortex knot
                        field_Ssb[id] = 0.f;
                }
                else {
                    // Far field value
                    field_Ssb[id] = 0.f;
                }

                if (field_Ssb[id] > 0.f) {
                    mean_sb_pos += field_Ssb[id];
                    Npos++;
                }
                else if (field_Ssb[id] < 0.f) {
                    mean_sb_neg += field_Ssb[id];
                    Nneg++;
                }


                if (field_Ssb[id] > max_sb) {
                    max_sb = field_Ssb[id];
                }
                if (field_Ssb[id] < min_sb) {
                    min_sb = field_Ssb[id];
                }
            }
        }
    }

    if (Npos > 0)
        mean_sb_pos /= (float)Npos;
    if (Nneg > 0)
        mean_sb_neg /= (float)Nneg;

    // Choose the cutoff values
    float up = 0.f;
    float un = 0.f;
    float alpha = 0.1f;
    float tc_pos = (1.f - up) * (mean_sb_pos/max_sb) * alpha + up;
    float tc_neg = (1.f - un) * (mean_sb_neg/min_sb) * alpha + un;


    // ===============================
    // Splay-bend pieces             |
    // ===============================
    // Positive splay-bend knot isosurface
    auto sub2ind_int = [&](int i, int j, int k) {
        return i + vox_interp[0] * j + slice_interp * k;
    };
    std::vector<uint> pos_splay_ribbon;
    for (const auto & id : fat_vortex_list) {

        int k = id / slice_interp;
        int j = (id - k * slice_interp) / vox_interp[0];
        int i = id - j * vox_interp[0] - k * slice_interp;
        //if (field_Ssb[id] > max_sb * _widget.knot_interaction_handle.sb_pos_iso_ratio_lower &&
        //    field_Ssb[id] <= max_sb * _widget.knot_interaction_handle.sb_pos_iso_ratio_upper) {
        //    pos_splay_ribbon.emplace_back(id);
        //}

        float local_avg_field_Ssb = field_Ssb[id];
        for (int x : {-1, 1})
            local_avg_field_Ssb += field_Ssb[sub2ind_int(i + x, j, k)];
        for (int y : {-1, 1})
            local_avg_field_Ssb += field_Ssb[sub2ind_int(i, j + y, k)];
        for (int z : {-1, 1})
            local_avg_field_Ssb += field_Ssb[sub2ind_int(i, j, k + z)];

        local_avg_field_Ssb /= 7.;

        //if (_widget.knot_interaction_handle.sb_pos_iso_ratio_lower * stdev < field_Ssb[id]
        //    && field_Ssb[id] < _widget.knot_interaction_handle.sb_pos_iso_ratio_upper * stdev) {
        //    pos_splay_ribbon.emplace_back(id);
        //}
                
        if (tc_pos <= local_avg_field_Ssb / max_sb)
            pos_splay_ribbon.emplace_back(id);

        //valid_fieldf[id] = 0.f;
    }
   
    // Reset sample grid
    for (auto id = 0; id < vol_interp; id++)
        sample_grid[id] = 0.f;

    // Choose the points to include in the isosurface
    for (const auto& id : pos_splay_ribbon) {
        sample_grid[id] = -1.f;
        //valid_fieldf[id] = -1.f;
    }


    mc2.reset_mesh();
    mc2.run();

    for (uint i = 0; i < mc2.nverts(); ++i) {
        LC::ExtendedMC::Vertex& v = mc2.vertices()[i];
        v.x = dr_interp[0] * v.x - cellf[0] * 0.5f;
        v.y = dr_interp[1] * v.y - cellf[1] * 0.5f;
        v.z = dr_interp[2] * v.z - cellf[2] * 0.5f;

        // Normalize normals
        float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
        if (nrm != 0)
        {
            nrm = 1.0 / sqrt(nrm);
            v.nx *= nrm;
            v.ny *= nrm;
            v.nz *= nrm;
        }
    }

    if (mc2.nverts() && mc2.ntrigs()) {
        unsigned int nVert = mc2.nverts();
        unsigned int nInd = mc2.ntrigs() * 3;
        LC::Math::IsoVertex* verts = new LC::Math::IsoVertex[nVert];
        unsigned int* indices = new unsigned int[nInd];

        LC::ExtendedMC::Vertex* temp_verts = mc2.vertices();
        int* ind_temp = (int*)mc2.triangles();

        // Feed data to indices,verts

        for (int i = 0; i < nInd; i++)
            indices[i] = ind_temp[i];

        for (int i = 0; i < nVert; i++) {
            verts[i].position[0] = temp_verts[i].x;
            verts[i].position[1] = temp_verts[i].y;
            verts[i].position[2] = temp_verts[i].z;
            verts[i].normal[0] = temp_verts[i].nx;
            verts[i].normal[1] = temp_verts[i].ny;
            verts[i].normal[2] = temp_verts[i].nz;

        }

        std::array<float, 4> blue = { 0.f,0.f,1.f,1.f };

        // Re-interleave data to IsoVertex format
        std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc2.nverts()]);
        LC::ExtendedMC::Vertex* iso_v_ptr = mc2.vertices();
        uint* iso_ind_ptr = (uint*)mc2.triangles();

        for (int vi = 0; vi < mc2.nverts(); vi++) {
            iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
        }


        // Smooth the mesh and save it
        //TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(), data->voxels, blue, false, 0, 0, saveObj, obj_name + std::string("_SB_pos.ply"));

        // Apply Taubin smoothing to the vertices
        float lambda = 0.35f;
        float mu = -0.34f;
        int smoothIterations = 80;
        auto smoothing_result = TaubinSmoothingGeneral(lambda, mu, smoothIterations, iso_verts.get(), mc2.nverts(), iso_ind_ptr,
            mc2.ntrigs() * 3, field_nn.get(), data->voxels, blue, false, 0, 0);


        smooth0_vertices = std::get<0>(smoothing_result);
        smooth0_trigs = std::get<1>(smoothing_result);

        Color4 c4_blue = { 0.f, 0.f, 1.f, 1.f }; // blue
        _vortexKnot.push_back(VortexKnot{});
        VortexKnot* knot;
        {
            // Go to the newly added vortex knot
            std::list<VortexKnot>::iterator iter;
            for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++)
                iter = it;
            knot = &(*iter);
        }

        knot->knotColor = c4_blue;
        knot->ribbon = true;

        // Load the surface to view
        // Initialize the mesh
        knot->surface.Init((LC::Surface::Vertex*)&smooth0_vertices[0], smooth0_vertices.size(), &smooth0_trigs[0][0],
            smooth0_trigs.size() * 3, _widget.preimage_translate);

        // Create the mesh
        knot->mesh = knot->surface.Mesh();

        // Add default vortex line to the scene
        new LC::Drawable::TransparentNormalDrawable{ *_vortexManipulator, _phongShader, *(knot->mesh), knot->draw, knot->cullFaces, _transparentNormalDrawables };
        mc2.reset_mesh();
    }
    else {
        LC_INFO("Not enough points for positive SB");
    }

    // yesn't
    // Negative splay bend
    std::vector<uint> neg_splay_ribbon;
    for (const auto& id : fat_vortex_list) {

        uint k = id / slice_interp;
        uint j = (id - k * slice_interp) / vox_interp[0];
        uint i = id - j * vox_interp[0] - k * slice_interp;
        // upper(mag) <= Ssb < lower(mag)
        //if (_widget.knot_interaction_handle.sb_neg_iso_ratio_upper * min_sb <= field_Ssb[id] &&
        //    field_Ssb[id] < min_sb * _widget.knot_interaction_handle.sb_neg_iso_ratio_lower) {
        //    neg_splay_ribbon.emplace_back(id);
        //}

        float local_avg_field_Ssb = field_Ssb[id];
        for (int x : {-1, 1})
            local_avg_field_Ssb += field_Ssb[sub2ind_int(i + x, j, k)];
        for (int y : {-1, 1})
            local_avg_field_Ssb += field_Ssb[sub2ind_int(i, j + y, k)];
        for (int z : {-1, 1})
            local_avg_field_Ssb += field_Ssb[sub2ind_int(i, j, k + z)];

        local_avg_field_Ssb /= 7.;


        //if (-_widget.knot_interaction_handle.sb_neg_iso_ratio_upper * stdev < field_Ssb[id]
        //    && field_Ssb[id] < -_widget.knot_interaction_handle.sb_neg_iso_ratio_lower * stdev) {
        //    neg_splay_ribbon.emplace_back(id);
        //}

        if (tc_neg <= local_avg_field_Ssb / min_sb)
            neg_splay_ribbon.emplace_back(id);

        //valid_fieldf[id] = 0.f;
    }

    // Reset sample grid
    for (auto id = 0; id < vol_interp; id++)
        sample_grid[id] = 0.f;

    // Choose the points to include in the isosurface
    for (const auto& id : neg_splay_ribbon) {
        sample_grid[id] = -1.f;
        //valid_fieldf[id] = -1.f;
    }

    mc2.reset_mesh();
    mc2.run();

    for (uint i = 0; i < mc2.nverts(); ++i) {
        LC::ExtendedMC::Vertex& v = mc2.vertices()[i];
        v.x = dr_interp[0] * v.x - cellf[0] * 0.5f;
        v.y = dr_interp[1] * v.y - cellf[1] * 0.5f;
        v.z = dr_interp[2] * v.z - cellf[2] * 0.5f;

        // Normalize normals
        float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
        if (nrm != 0)
        {
            nrm = 1.0 / sqrt(nrm);
            v.nx *= nrm;
            v.ny *= nrm;
            v.nz *= nrm;
        }
    }

    if (mc2.nverts() && mc2.ntrigs()) {
        unsigned int nVert = mc2.nverts();
        unsigned int nInd = mc2.ntrigs() * 3;
        LC::Math::IsoVertex* verts = new LC::Math::IsoVertex[nVert];
        unsigned int* indices = new unsigned int[nInd];

        LC::ExtendedMC::Vertex* temp_verts = mc2.vertices();
        int* ind_temp = (int*)mc2.triangles();

        // Feed data to indices,verts

        for (int i = 0; i < nInd; i++)
            indices[i] = ind_temp[i];

        for (int i = 0; i < nVert; i++) {
            verts[i].position[0] = temp_verts[i].x;
            verts[i].position[1] = temp_verts[i].y;
            verts[i].position[2] = temp_verts[i].z;
            verts[i].normal[0] = temp_verts[i].nx;
            verts[i].normal[1] = temp_verts[i].ny;
            verts[i].normal[2] = temp_verts[i].nz;

        }


        std::vector<MeshLib::PNCVertex<float>> vertices;
        std::vector<MeshLib::Triangle> triangles;
        std::array<float, 4> yellow = { 1.f,1.f,0.f,1.f };

        // Re-interleave data to IsoVertex format
        std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc2.nverts()]);
        LC::ExtendedMC::Vertex* iso_v_ptr = mc2.vertices();
        uint* iso_ind_ptr = (uint*)mc2.triangles();

        for (int vi = 0; vi < mc2.nverts(); vi++) {
            iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
        }

        // Smooth the mesh and save it
        //TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(), data->voxels, yellow, false, 0, 0, saveObj, obj_name + std::string("_SB_neg.ply"));

        float lambda = 0.33f;
        float mu = -0.34f;
        int smoothIterations = 80;
        auto smoothing_result = TaubinSmoothingGeneral(lambda, mu, smoothIterations, iso_verts.get(), mc2.nverts(), iso_ind_ptr,
            mc2.ntrigs() * 3, field_nn.get(), data->voxels, yellow, false, 0, 0);

        smooth0_vertices = std::get<0>(smoothing_result);
        smooth0_trigs = std::get<1>(smoothing_result);

        Color4 c4_yellow = { 1.f, 1.f, 0.f, 1.f }; // blue
        _vortexKnot.push_back(VortexKnot{});
        VortexKnot* knot;
        {
            // Go to the newly added vortex knot
            std::list<VortexKnot>::iterator iter;
            for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++)
                iter = it;
            knot = &(*iter);
        }

        knot->knotColor = c4_yellow;
        knot->ribbon = true;

        // Load the surface to view
        // Initialize the mesh
        knot->surface.Init((LC::Surface::Vertex*)&smooth0_vertices[0], smooth0_vertices.size(), &smooth0_trigs[0][0],
            smooth0_trigs.size() * 3, _widget.preimage_translate);

        // Create the mesh
        knot->mesh = knot->surface.Mesh();

        // Add default vortex line to the scene
        new LC::Drawable::TransparentNormalDrawable{ *_vortexManipulator, _phongShader, *(knot->mesh), knot->draw, knot->cullFaces, _transparentNormalDrawables };
        mc2.reset_mesh();
    }
    else {
        LC_INFO("Not enough points for negative SB");
    }


}

std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>>
Sandbox::findVortexKnot(bool saveObj, const std::string& obj_name,
    // Parameters
    std::unique_ptr<float[]> &chi_field, std::unique_ptr<float[]>& field_nn, std::unique_ptr<float[]>& valid_fieldf,
    std::unique_ptr<float[]>& field_S, std::unique_ptr<float[]>& sample_grid,
    std::unique_ptr<float[]>& field_Q, std::unique_ptr<float[]>& field_Ssb,
    LC::ExtendedMC::MarchingCubes &mc,
    LC::ExtendedMC::MarchingCubes &mc2,
    std::array<int,3> &vox_interp
    ) {
    Dataset* data = (Dataset*)(_solver->GetDataPtr());
    std::array<float, 3> cellf = { data->cell_dims[0], data->cell_dims[1], data->cell_dims[2] };
    std::array<LC::scalar, 3> dr;
    for (int i = 0; i < 3; i++)
        dr[i] = cellf[i] / (data->voxels[i] - 1);


    // Parameters for data management
    uint slice = data->voxels[0] * data->voxels[1];
    uint vol = data->voxels[2] * slice;
    std::array<float, 3> dr_interp;
    for (int d = 0; d < 3; d++) {
        dr_interp[d] = cellf[d] / (vox_interp[d] - 1);
    }
    uint vol_interp = vox_interp[0] * vox_interp[1] * vox_interp[2];

    // Clean knots
    if (!_vortexKnot.empty()) {
        VortexKnot* knot;
        for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++) {
            knot = &(*it);
            knot->draw = false;
        }

        while (!_vortexKnot.empty()) {
            std::list<VortexKnot>::iterator iter;
            for (auto it = _vortexKnot.begin(); it != _vortexKnot.end(); it++)
                iter = it;
            _vortexKnot.erase(iter);
        }
    }


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

    // Compute everything
    LC::Math::ChiralityField(field_nn.get(), chi_field, data->voxels, cellf, _vortex_line.valid_field, false, _widget.knot_interaction_handle.tolerance, false);
    

    // Convert valid field to floats and make a list of all defect indices
    std::vector<uint> vortex_indices;
    std::set<uint> fat_vortex_list;

    for (auto ii = 0; ii < vol; ii++) {

        valid_fieldf[ii] = 0.f;
        if (_vortex_line.valid_field[ii]) {

            int k = ii / slice;
            int j = (ii - k * slice) / data->voxels[0];
            int i = ii - j * data->voxels[1] - slice * k;
            LC::scalar avgcos2_chi = 0.;
            LC::scalar avgcos2_tau = 0.;
            Eigen::Vector3d chi0(chi_field[ii], chi_field[ii + vol], chi_field[ii + 2 * vol]);
            Eigen::Vector3d nn0(field_nn[ii], field_nn[ii + vol], field_nn[ii + 2 * vol]);
            Eigen::Vector3d tau0 = chi0.cross(nn0);
            if (tau0.squaredNorm() > 0.)
                tau0.normalize();
            // Compute the scalar order parameter with the chi field


            int points_in_bds = 0;

            for (int x = -1; x <= 1; x++) {
                for (int y = -1; y <= 1; y++) {
                    for (int z = -1; z <= 1; z++) {
                        uint ip = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                        if (ip >= 0 && ip < vol && ip != ii) {
                            points_in_bds++;
                            Eigen::Vector3d chi_xyz(chi_field[ip], chi_field[ip + vol], chi_field[ip + 2 * vol]);
                            Eigen::Vector3d nn_xyz(field_nn[ip], field_nn[ip + vol], field_nn[ip + 2 * vol]);
                            Eigen::Vector3d tau_xyz = chi_xyz.cross(nn_xyz);

                            if (tau_xyz.squaredNorm() > 0.)
                                tau_xyz.normalize();

                            avgcos2_chi += pow(chi0.dot(chi_xyz), 2);
                            avgcos2_tau += pow(tau0.dot(tau_xyz), 2);
                            
                        }
                    }
                }
            }

            LC::scalar Sop_chi, Sop_lambda, Sop_tau;
            if (points_in_bds > 0) {
                Sop_chi = 0.5 * (3. * avgcos2_chi / points_in_bds - 1.);
                Sop_tau = 0.5 * (3. * avgcos2_tau / points_in_bds - 1.);
                LC::scalar S = 0.5 * (Sop_chi + Sop_tau);
                field_S[ii] = S;
                if (abs(S) < _widget.knot_interaction_handle.isoValue) {
                    //vortex_indices.emplace_back(ii);
                    _vortex_line.valid_field[ii] = 0;
                }
            }

        }
        else {
            field_S[ii] = 0.;
            //vortex_indices.emplace_back(ii);
        }
    }

    // New vortex classification scheme
    for (auto ii = 0; ii < vol; ii++) {
        int k = ii / slice;
        int j = (ii - k * slice) / data->voxels[0];
        int i = ii - j * data->voxels[1] - slice * k;

        int line = data->voxels[0];

        auto get_id = [slice, line](int x, int y, int z) {
            return x + line * y + slice * z;
        };

        // Compute S derivatives
        float Sd2 = 0.f;
        if (i > 1 && j > 1 && k > 1 && i < data->voxels[0] - 2 && j < data->voxels[1] - 2 && k < data->voxels[2] - 2) {
            float S200 = (field_S[get_id(i + 1, j, k)] + field_S[get_id(i - 1, j, k)] - 2. * field_S[ii]);
            float S020 = (field_S[get_id(i, j + 1, k)] + field_S[get_id(i, j - 1, k)] - 2. * field_S[ii]);
            float S002 = (field_S[get_id(i, j, k + 1)] + field_S[get_id(i, j, k - 1)] - 2. * field_S[ii]);
            Sd2 = S200 + S020 + S002;
        }


        // Classify the current index ii
        if (field_S[ii] < _widget.knot_interaction_handle.isoValue) {
            //if (Sd2 > 0.5)
                vortex_indices.emplace_back(ii);
        }

    }

    // Classify only the raw defect core
    for (const auto& id : vortex_indices) {
        valid_fieldf[id] = -1.f;
    }

    // Recount newly classified nodes
    _vortex_line.numNodes = 0;

    for (auto ii = 0; ii < vol; ii++) {
        if (valid_fieldf[ii] < 0.f) {
            _vortex_line.valid_field[ii] = 0;
            _vortex_line.numNodes++;
        }
    }

    std::map<uint, Face> unvisited_list_local;
    // Key: vertex, value: Triangle
    std::multimap<uint, uint> mmap;
    // Make the initial color white to distinguish components that were too small
    std::array<float, 4> col_white = { 1.f, 1.f, 1.f, 1.f };

    auto extended_to_iso = [](const LC::ExtendedMC::Vertex& v) {
        LC::Math::IsoVertex viso;
        viso.position[0] = v.x; viso.position[1] = v.y; viso.position[2] = v.z;
        viso.normal[0] = v.nx; viso.normal[1] = v.ny; viso.normal[2] = v.nz;

        return viso;
    };

    // Use extended marching cubes to generate the isosurface

    mc.set_ext_data(valid_fieldf.get());
    mc2.set_ext_data(sample_grid.get());
    mc.reset_mesh();
    mc2.reset_mesh();

    mc.run();

    std::tuple<std::vector<Eigen::Vector3d>, std::vector<int>> exported_data; // exported vortex knot data

    // Smooth the raw vortex knot isosurface
    std::vector<MeshLib::PNCVertex<float>> smooth0_vertices;
    std::vector<MeshLib::Triangle> smooth0_trigs;
    {
        // Re-interleave data to IsoVertex format
        std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc.nverts()]);
        LC::ExtendedMC::Vertex* iso_v_ptr = mc.vertices();
        uint* iso_ind_ptr = (uint*)mc.triangles();

        for (int vi = 0; vi < mc.nverts(); vi++) {
            iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
        }

        LC::ExtendedMC::Triangle* tri_data = mc.triangles();

        // Random color
        std::array<float, 4> col_arr;
        for (int d = 0; d < 4; d++)
            col_arr[d] = col_white[d];

        // Apply Taubin smoothing to the pre-sampled vertices
        float lambda = 0.5f;
        float mu = -0.34f;
        int smoothIterations = 80;

        auto smoothing_result = TaubinSmoothingGeneral(lambda, mu, smoothIterations, iso_verts.get(), mc.nverts(), iso_ind_ptr,
            mc.ntrigs() * 3, 0, data->voxels, col_arr, 0, 0, 0, 0);

        smooth0_vertices = std::get<0>(smoothing_result);
        smooth0_trigs = std::get<1>(smoothing_result);

        // Clean the mc data
        mc.reset_mesh();
    }

    // Rescale positions
    {
        float dx = data->cell_dims[0] / (data->voxels[0] - 1);
        float dy = data->cell_dims[1] / (data->voxels[1] - 1);
        float dz = data->cell_dims[2] / (data->voxels[2] - 1);
        float xmin = -0.5 * data->cell_dims[0];
        float ymin = -0.5 * data->cell_dims[1];
        float zmin = -0.5 * data->cell_dims[2];

        for (uint i = 0; i < smooth0_vertices.size(); ++i) {
            MeshLib::PNCVertex<float>& v = smooth0_vertices[i];

            v.position[0] = dx * v.position[0] + xmin;
            v.position[1] = dy * v.position[1] + ymin;
            v.position[2] = dz * v.position[2] + zmin;

            // Normalize normals
            float nrm = v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] + v.normal[2] * v.normal[2];
            if (nrm != 0)
            {
                nrm = 1.0 / sqrt(nrm);
                v.normal[0] *= nrm;
                v.normal[1] *= nrm;
                v.normal[2] *= nrm;
            }
        }
    }
    

    if (smooth0_vertices.size() && smooth0_trigs.size()) {
        unsigned int nVert = smooth0_vertices.size();
        unsigned int nInd = smooth0_trigs.size() * 3;
        unsigned int nTriangles = nInd / 3;

        // Feed data to indices,verts

        //mc.reset_mesh();

        // Note that verts currently refers to verts and not global indices
        // Therefore, convert indices data to global index space through vertex positions
        for (int tri = 0; tri < nTriangles; tri++) {
            uint i1 = smooth0_trigs[tri][0];
            uint i2 = smooth0_trigs[tri][1];
            uint i3 = smooth0_trigs[tri][2];

            mmap.insert({ i1, tri });
            mmap.insert({ i2, tri });
            mmap.insert({ i3, tri });

            unvisited_list_local.insert({ tri, Face(i1, i2, i3) });
        }

        // Find components through the vertices
        std::vector<std::vector<uint>> components;
        try {
            components = find_all_components_graph(unvisited_list_local, nVert, 15);
        }
        catch (const std::runtime_error& e) {
            LC_INFO("Error: {0}", e.what());
        }
        
        std::array<float, 4> randcol = { 1.f,1.f,1.f,1.f };
        std::tuple <std::vector<MeshLib::PNCVertex<float>>, std::vector<MeshLib::Triangle>> smoothing_result;
        
        // Determine the minimum point separation of vortex knot points
        LC::scalar dist = 2. * _widget.knot_interaction_handle.point_density * sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

        // Export string knots
        if (_widget.knot_interaction_handle.saveStringKnots) {
            exported_data = exportVortexKnot(smooth0_vertices, smooth0_trigs, dist, obj_name + std::string("_vortex.bin"));
            // Transform positions to director data
            auto pdata = std::get<0>(exported_data);
            std::vector < Eigen::Vector2d > odata(pdata.size());

            for (int p = 0; p < pdata.size(); p++) {
                // Transform position to index coords
                int i = (pdata[p].x() / cellf[0] + 0.5) * (data->voxels[0] - 1);
                int j = (pdata[p].y() / cellf[1] + 0.5) * (data->voxels[1] - 1);
                int k = (pdata[p].z() / cellf[2] + 0.5) * (data->voxels[2] - 1);
                uint id = i + j * data->voxels[0] + k * data->voxels[0] * data->voxels[1];

                odata[p][0] = acos(field_nn[id + 2 * vol]); // theta
                odata[p][1] = atan2(field_nn[id + vol], field_nn[id]); // phi
            }

            // Export
            std::string fname = obj_name + "_sphere.bin";
            std::ofstream ofile(fname.c_str(), std::ios::out | std::ios::binary);
            
            if (ofile.is_open()) {
                int32_t osz = odata.size();
                ofile.write((char*)&osz, sizeof(int32_t));
                ofile.write((char*)&odata[0], osz * sizeof(Eigen::Vector2d));
                ofile.close();
            }
        }

        int Radius = _widget.knot_interaction_handle.saturation;

        // Split apart vertex data
        int comp = 0;
        for (auto& component : components) {

            // Initialize the sample grid
            for (int ii = 0; ii < vol_interp; ii++)
                sample_grid[ii] = 0.1f;

            // Sample the component using the indices
            for (const auto& ci : component) {
                auto v_ci = smooth0_vertices[ci];

                // Convert the vertex into an index on the new grid
                uint i = (v_ci.position[0] + 0.5f * cellf[0]) / cellf[0] * (vox_interp[0] - 1);
                uint j = (v_ci.position[1] + 0.5f * cellf[1]) / cellf[1] * (vox_interp[1] - 1);
                uint k = (v_ci.position[2] + 0.5f * cellf[2]) / cellf[2] * (vox_interp[2] - 1);

                // Calculate the new index
                uint id = i + j * vox_interp[0] + k * vox_interp[0] * vox_interp[1];

                // Apply some sampling technique to the vortex knot at this position
                for (int x = -Radius; x <= Radius; x++)
                    for (int y = -Radius; y <= Radius; y++)
                        for (int z = -Radius; z <= Radius; z++) {
                            // Compute the index to modify
                            uint id_ball = (i + x) + (j + y) * vox_interp[0] + (k + z) * vox_interp[0] * vox_interp[1];
                            // Make sure defect is in bounds
                            if (id_ball < vol_interp && id_ball >= 0) {
                                //_vortex_line.valid_field[id_ball] = 0;
                                float r2 = x * x + y * y + z * z;
                                if (id_ball == id)
                                    sample_grid[id_ball] += -1.f;
                                else
                                    sample_grid[id_ball] += -1.f / r2;
                            }
                        }

            }

            // Compute the new knot isosurface from sample_grid

            mc2.reset_mesh();
            mc2.run();

            // Rescale the new isosurface vertices
            {
                float xmin = -0.5 * cellf[0];
                float ymin = -0.5 * cellf[1];
                float zmin = -0.5 * cellf[2];

                for (uint i = 0; i < mc2.nverts(); ++i) {
                    LC::ExtendedMC::Vertex& v = mc2.vertices()[i];
                    v.x = dr_interp[0] * v.x + xmin;
                    v.y = dr_interp[1] * v.y + ymin;
                    v.z = dr_interp[2] * v.z + zmin;

                    // Normalize normals
                    float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
                    if (nrm != 0)
                    {
                        nrm = 1.0 / sqrt(nrm);
                        v.nx *= nrm;
                        v.ny *= nrm;
                        v.nz *= nrm;
                    }
                }
            }

            // Re-interleave data to IsoVertex format
            std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc2.nverts()]);
            LC::ExtendedMC::Vertex* iso_v_ptr = mc2.vertices();
            uint* iso_ind_ptr = (uint*)mc2.triangles();

            for (int vi = 0; vi < mc2.nverts(); vi++) {
                iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
            }

            LC::ExtendedMC::Triangle* tri_data = mc2.triangles();

            // Apply Taubin smoothing to the vertices
            float lambda = 0.5f;
            float mu = -0.34f;
            int smoothIterations = 120;
            smoothing_result = TaubinSmoothingGeneral(lambda,mu,smoothIterations, iso_verts.get(), mc2.nverts(), iso_ind_ptr,
                mc2.ntrigs() * 3, field_nn.get(), data->voxels, randcol, false, 1, 0, true, 
                obj_name + std::string("_vortex"+std::to_string(++comp)+"_of_"+ std::to_string(components.size()) +".ply"));

            
            
            //auto smooth_vertices = std::get<0>(smoothing_result);
            //auto smooth_triangles = std::get<1>(smoothing_result);
        }
    }

    //std::tuple<std::vector<Eigen::Vector3d>, std::vector<int>> exported_data; // exported vortex knot data
    //exported_data = exportVortexKnot(vertices, triangles, dist, obj_name + std::string("_vortex.bin"));
    //smoothing_result = TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(), data->voxels, randcol, false, 1, 0, true, obj_name + std::string("_vortex.ply"));
#define SAVE_SPLAY_BEND 1

#if SAVE_SPLAY_BEND
    int Nx = data->voxels[0];

    auto sub2ind = [slice, Nx, vol](int i, int j, int k, int d) {
        return i + Nx * j + slice * k + vol * d;
    };
    // Compute Qij for the chi field
    for (int i = 0; i < data->voxels[0]; i++) {
        for (int j = 0; j < data->voxels[1]; j++) {
            for (int k = 0; k < data->voxels[2]; k++) {

                uint ii = i + data->voxels[0] * j + slice * k;
                // Compute the order parameter S

                Eigen::Vector3d chi0(chi_field[ii], chi_field[ii + vol], chi_field[ii + 2 * vol]);
                // Compute the scalar order parameter with the chi field

                // Use the scalar order parameter value computed above
                LC::scalar S = 1.;// field_S[ii];


                field_Q[ii] = 0.5 * S * (3. * chi0[0] * chi0[0] - 1.); // qa
                field_Q[ii + vol] = 0.5 * S * (3. * chi0[0] * chi0[1]); // qb
                field_Q[ii + 2 * vol] = 0.5 * S * (3. * chi0[0] * chi0[2]); // qc
                field_Q[ii + 3 * vol] = 0.5 * S * (3. * chi0[1] * chi0[1] - 1.); // qd
                field_Q[ii + 4 * vol] = 0.5 * S * (3. * chi0[1] * chi0[2]); // qe
            }
        }
    }

    // Fill fattened vortex
    int sb_rad = _widget.knot_interaction_handle.sb_saturation;
    for (const auto& id : vortex_indices) {
        // Convert the index to i,j,k space
        int k = id / slice;
        int j = (id - k * slice) / data->voxels[0];
        int i = id - j * data->voxels[0] - slice * k;

        for (int x = -sb_rad; x <= sb_rad; x++)
            for (int y = -sb_rad; y <= sb_rad; y++)
                for (int z = -sb_rad; z <= sb_rad; z++) {
                    // Compute the index to modify
                    uint id_ball = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                    // Make sure defect is in bounds
                    if (id_ball < vol && id_ball >= 0) {
                        //_vortex_line.valid_field[id_ball] = 0;
                        fat_vortex_list.insert(id_ball);
                    }
                }

    }

    // Statistical parameters
    float max_sb = 0.f;
    float min_sb = 0.f;
    float mean = 0.f;
    float stdev = 0.f;
    int N = 0;

    // [derivative][qtensor]
    //float qD[3][5];
    // [position][qtensor][--, -, +, ++]
    //float Qt[3][5][4];
    //constexpr float c1 = 1.0 / 12.0;
    //constexpr float c2 = 2.0 / 3.0;
    //constexpr float c3 = 4.0 / 3.0;

    //auto D2 = [&](int i, int d) {
    //    qD[i][d] = (-c1 * Qt[i][d][3] + c2 * Qt[i][d][2] - c2 * Qt[i][d][1] + c1 * Qt[i][d][0]) / dr[i];
    //};

    // Next compute the seven terms needed for splay bend
    for (int i = 0; i < data->voxels[0]; i++) {
        for (int j = 0; j < data->voxels[1]; j++) {
            for (int k = 0; k < data->voxels[2]; k++) {
                if (i > 1 && i < data->voxels[0] - 2 && j > 1 && j < data->voxels[1] - 2 && k > 1 && k < data->voxels[2] - 2) {
                    
                    // Second order
                    float qa002 = (field_Q[sub2ind(i, j, k + 1, 0)] + field_Q[sub2ind(i, j, k - 1, 0)] - 2.f * field_Q[sub2ind(i, j, k, 0)]) / (dr[2] * dr[2]);
                    float qa200 = (field_Q[sub2ind(i + 1, j, k, 0)] + field_Q[sub2ind(i - 1, j, k, 0)] - 2.f * field_Q[sub2ind(i, j, k, 0)]) / (dr[0] * dr[0]);
                    float qd002 = (field_Q[sub2ind(i, j, k + 1, 3)] + field_Q[sub2ind(i, j, k - 1, 3)] - 2.f * field_Q[sub2ind(i, j, k, 3)]) / (dr[2] * dr[2]);
                    float qd020 = (field_Q[sub2ind(i, j + 1, k, 3)] + field_Q[sub2ind(i, j - 1, k, 3)] - 2.f * field_Q[sub2ind(i, j, k, 3)]) / (dr[1] * dr[1]);
                    float qb110 = (field_Q[sub2ind(i + 1, j + 1, k, 1)] + field_Q[sub2ind(i - 1, j - 1, k, 1)]
                        - field_Q[sub2ind(i - 1, j + 1, k, 1)] - field_Q[sub2ind(i + 1, j - 1, k, 1)]) / (4.f * dr[0] * dr[1]);
                    float qc101 = (field_Q[sub2ind(i + 1, j, k + 1, 2)] + field_Q[sub2ind(i - 1, j, k - 1, 2)]
                        - field_Q[sub2ind(i - 1, j, k + 1, 2)] - field_Q[sub2ind(i + 1, j, k - 1, 2)]) / (4.f * dr[0] * dr[2]);
                    float qe011 = (field_Q[sub2ind(i, j + 1, k + 1, 4)] + field_Q[sub2ind(i, j - 1, k - 1, 4)]
                        - field_Q[sub2ind(i, j - 1, k + 1, 4)] - field_Q[sub2ind(i, j + 1, k - 1, 4)]) / (4.f * dr[1] * dr[2]);
                    

// Doesn't work very well because the derivatives are highly localized
#if 0

                    for (int d = 0; d < 5; d++) {

                        // Fill
                        Qt[0][d][0] = field_Q[sub2ind(i - 2, j, k, d)];
                        Qt[0][d][1] = field_Q[sub2ind(i - 1, j, k, d)];
                        Qt[0][d][2] = field_Q[sub2ind(i + 1, j, k, d)];
                        Qt[0][d][3] = field_Q[sub2ind(i + 2, j, k, d)];

                        Qt[1][d][0] = field_Q[sub2ind(i, j - 2, k, d)];
                        Qt[1][d][1] = field_Q[sub2ind(i, j - 1, k, d)];
                        Qt[1][d][2] = field_Q[sub2ind(i, j + 1, k, d)];
                        Qt[1][d][3] = field_Q[sub2ind(i, j + 2, k, d)];

                        Qt[2][d][0] = field_Q[sub2ind(i, j, k - 2, d)];
                        Qt[2][d][1] = field_Q[sub2ind(i, j, k - 1, d)];
                        Qt[2][d][2] = field_Q[sub2ind(i, j, k + 1, d)];
                        Qt[2][d][3] = field_Q[sub2ind(i, j, k + 2, d)];

                        for (int i = 0; i < 3; i++) {
                            qD[i][d] = (-c1 * Qt[i][d][3] + c2 * Qt[i][d][2] - c2 * Qt[i][d][1] + c1 * Qt[i][d][0]) / dr[i];
                        }
                        D2(2, 0);
                        D2(0, 0);
                        D2(2, 3);
                        D2(1, 3);
                    }
                    // fourth order
                    float qa002 = qD[2][0];
                    float qa200 = qD[0][0];
                    float qd002 = qD[2][3];
                    float qd020 = qD[1][3];
                        // (-2, -2) (2, 2) (-2, 2) (2, -2)
                        //    +       +       -       -
                    float qb110 = (field_Q[sub2ind(i-2, j-2, k,1)] + field_Q[sub2ind(i+2, j+2, k,1)] - field_Q[sub2ind(i-2, j+2, k,1)] - field_Q[sub2ind(i+2, j-2, k,1)] +

                        // (1, -2) (-2, 1) (-1, -2) (-2, -1)
                        //    +       +       -       -
                        (field_Q[sub2ind(i+1, j-2, k,1)] + field_Q[sub2ind(i-2, j+1, k,1)] - field_Q[sub2ind(i-1, j-2, k,1)] - field_Q[sub2ind(i-2, j-1, k,1)] +

                        // (2, -1) (-1, 2) (1, 2) (2, 1)
                        //    +       +       -       -
                        field_Q[sub2ind(i+2, j-1, k,1)] + field_Q[sub2ind(i-1, j+2, k,1)] - field_Q[sub2ind(i+1, j+2, k,1)] - field_Q[sub2ind(i+2, j+1, k,1)]) * 8.0f +

                        // (1, 1) (-1, -1) (1, -1) (-1, 1)
                        //    +       +       -       -
                        (field_Q[sub2ind(i+1, j+1, k,1)] + field_Q[sub2ind(i-1, j-1, k,1)] - field_Q[sub2ind(i+1, j-1, k,1)] - field_Q[sub2ind(i-1, j+1, k,1)]) * 64.0f) / (144.0f * dx * dy);


                    // (-2, -2) (2, 2) (-2, 2) (2, -2)
                        //    +       +       -       -
                    float qc101 = (field_Q[sub2ind(i-2, j, k-2,2)] + field_Q[sub2ind(i+2, j, k+2,2)] - field_Q[sub2ind(i-2, j, k+2,2)] - field_Q[sub2ind(i+2, j, k-2,2)] +

                        // (1, -2) (-2, 1) (-1, -2) (-2, -1)
                        //    +       +       -       -
                        (field_Q[sub2ind(i+1, j, k-2,2)] + field_Q[sub2ind(i-2, j, k+1,2)] - field_Q[sub2ind(i-1, j, k-2,2)] - field_Q[sub2ind(i-2, j, k-1,2)] +

                            // (2, -1) (-1, 2) (1, 2) (2, 1)
                            //    +       +       -       -
                            field_Q[sub2ind(i+2, j, k-1,2)] + field_Q[sub2ind(i-1, j, k+2,2)] - field_Q[sub2ind(i+1, j, k+2,2)] - field_Q[sub2ind(i+2, j, k+1,2)]) * 8.0f +

                        // (1, 1) (-1, -1) (1, -1) (-1, 1)
                        //    +       +       -       -
                        (field_Q[sub2ind(i+1, j, k+1,2)] + field_Q[sub2ind(i-1, j, k-1,2)] - field_Q[sub2ind(i+1, j, k-1,2)] - field_Q[sub2ind(i-1, j, k+1,2)]) * 64.0f) / (144.0f * dx * dz);



                    // (-2, -2) (2, 2) (-2, 2) (2, -2)
                    //    +       +       -       -

                    float qe011 = (field_Q[sub2ind(i, j-2, k-2,4)] + field_Q[sub2ind(i, j+2, k+2,4)] - field_Q[sub2ind(i, j-2, k+2,4)] - field_Q[sub2ind(i, j+2, k-2,4)] +

                        // (1, -2) (-2, 1) (-1, -2) (-2, -1)
                        //    +       +       -       -
                        (field_Q[sub2ind(i, j+1, k-2,4)] + field_Q[sub2ind(i, j-2, k+1,4)] - field_Q[sub2ind(i, j-1, k-2,4)] - field_Q[sub2ind(i, j-2, k-1,4)] +

                            // (2, -1) (-1, 2) (1, 2) (2, 1)
                            //    +       +       -       -
                            field_Q[sub2ind(i, j+2, k-1,4)] + field_Q[sub2ind(i, j-1, k+2,4)] - field_Q[sub2ind(i, j+1, k+2,4)] - field_Q[sub2ind(i, j+2, k+1,4)]) * 8.0f +

                        // (1, 1) (-1, -1) (1, -1) (-1, 1)          
                        //    +       +       -       -
                        (field_Q[sub2ind(i, j+1, k+1,4)] + field_Q[sub2ind(i, j-1, k-1,4)] - field_Q[sub2ind(i, j+1, k-1,4)] - field_Q[sub2ind(i, j-1, k+1,4)]) * 64.0f) / (144.0f * dy * dz);
#endif

                    uint id = sub2ind(i, j, k, 0);

                    // Found area to compute the splay bend
                    if (fat_vortex_list.find(id) != fat_vortex_list.end())
                        field_Ssb[id] = -qa002 + qa200 + 2.f * qb110 + 2.f * qc101 - qd002 + qd020 + 2.f * qe011;
                    else // splay bend not close enough to the vortex knot
                        field_Ssb[id] = 0.f;
                }
                else {
                    // Far field value
                    field_Ssb[sub2ind(i, j, k, 0)] = 0.f;
                }

           
                mean += field_Ssb[sub2ind(i, j, k, 0)];
                N++;
                

                if (field_Ssb[sub2ind(i, j, k, 0)] > max_sb) {
                    max_sb = field_Ssb[sub2ind(i, j, k, 0)];
                }
                if (field_Ssb[sub2ind(i, j, k, 0)] < min_sb) {
                    min_sb = field_Ssb[sub2ind(i, j, k, 0)];
                }
            }
        }
    }


    mean /= (float)N;


    // Compute the standard deviation
    for (auto i = 0; i < vol; i++) {
        stdev += pow(field_Ssb[i] - mean,2);
    }

    if (N > 0)
        stdev = sqrt(stdev / N);

    // Compute the two isosurfaces for positive and negative-splay bend
    // -----------------------------------------------------------------
    // > Positive splay bend
    std::vector<uint> pos_splay_ribbon;
    for (int i = 1; i < data->voxels[0]-1; i++) {
        for (int j = 1; j < data->voxels[1]-1; j++) {
            for (int k = 1; k < data->voxels[2]-1; k++) {
                uint id = sub2ind(i, j, k, 0);
                //if (field_Ssb[id] > max_sb * _widget.knot_interaction_handle.sb_pos_iso_ratio_lower &&
                //    field_Ssb[id] <= max_sb * _widget.knot_interaction_handle.sb_pos_iso_ratio_upper) {
                //    pos_splay_ribbon.emplace_back(id);
                //}

                float local_avg_field_Ssb = field_Ssb[id];
                for (int x : {-1, 1})
                    local_avg_field_Ssb += field_Ssb[sub2ind(i + x, j, k, 0)];
                for (int y : {-1, 1})
                    local_avg_field_Ssb += field_Ssb[sub2ind(i, j + y, k, 0)];
                for (int z : {-1, 1})
                    local_avg_field_Ssb += field_Ssb[sub2ind(i, j, k + z, 0)];

                local_avg_field_Ssb /= 7.;

                if (_widget.knot_interaction_handle.sb_pos_iso_ratio_lower* stdev < field_Ssb[id]
                    && field_Ssb[id] < _widget.knot_interaction_handle.sb_pos_iso_ratio_upper * stdev) {
                    pos_splay_ribbon.emplace_back(id);
                }
                valid_fieldf[id] = 0.f;
            }
        }
    }

    int Radius = 0;

    // Saturate the ribbon to enhance the contrast
    for (const auto& id : pos_splay_ribbon) {
        // Convert the index to i,j,k space
        int k = id / slice;
        int j = (id - k * slice) / data->voxels[0];
        int i = id - j * data->voxels[0] - slice * k;

        for (int x = -Radius; x <= Radius; x++)
            for (int y = -Radius; y <= Radius; y++)
                for (int z = -Radius; z <= Radius; z++) {
                    // Compute the index to modify
                    uint id_ball = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                    // Make sure defect is in bounds
                    if (id_ball < vol && id_ball >= 0) {
                        //_vortex_line.valid_field[id_ball] = 0;
                        float r2 = x * x + y * y + z * z;
                        valid_fieldf[id_ball] = -1.f;
                    }
                }

    }

    mc.reset_mesh();
    mc.run();

    for (uint i = 0; i < mc.nverts(); ++i) {
        LC::ExtendedMC::Vertex& v = mc.vertices()[i];
        v.x = dr[0] * v.x - cellf[0] * 0.5f;
        v.y = dr[1] * v.y - cellf[1] * 0.5f;
        v.z = dr[2] * v.z - cellf[2] * 0.5f;

        // Normalize normals
        float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
        if (nrm != 0)
        {
            nrm = 1.0 / sqrt(nrm);
            v.nx *= nrm;
            v.ny *= nrm;
            v.nz *= nrm;
        }
    }

    if (mc.nverts() && mc.ntrigs()) {
        unsigned int nVert = mc.nverts();
        unsigned int nInd = mc.ntrigs() * 3;
        LC::Math::IsoVertex* verts = new LC::Math::IsoVertex[nVert];
        unsigned int* indices = new unsigned int[nInd];

        LC::ExtendedMC::Vertex* temp_verts = mc.vertices();
        int* ind_temp = (int*)mc.triangles();

        // Feed data to indices,verts

        for (int i = 0; i < nInd; i++)
            indices[i] = ind_temp[i];

        for (int i = 0; i < nVert; i++) {
            verts[i].position[0] = temp_verts[i].x;
            verts[i].position[1] = temp_verts[i].y;
            verts[i].position[2] = temp_verts[i].z;
            verts[i].normal[0] = temp_verts[i].nx;
            verts[i].normal[1] = temp_verts[i].ny;
            verts[i].normal[2] = temp_verts[i].nz;

        }

        std::array<float, 4> blue = { 0.f,0.f,1.f,1.f };

        // Re-interleave data to IsoVertex format
        std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc.nverts()]);
        LC::ExtendedMC::Vertex* iso_v_ptr = mc.vertices();
        uint* iso_ind_ptr = (uint*)mc.triangles();

        for (int vi = 0; vi < mc.nverts(); vi++) {
            iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
        }


        // Smooth the mesh and save it
        //TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(), data->voxels, blue, false, 0, 0, saveObj, obj_name + std::string("_SB_pos.ply"));

        // Apply Taubin smoothing to the vertices
        float lambda = 0.35f;
        float mu = -0.34f;
        int smoothIterations = 80;
        TaubinSmoothingGeneral(lambda, mu, smoothIterations, iso_verts.get(), mc.nverts(), iso_ind_ptr,
            mc.ntrigs() * 3, field_nn.get(), data->voxels, blue, false, 0, 0, saveObj,
            obj_name + std::string("_SB_pos.ply"));

        mc.reset_mesh();
    }

    // > Negative splay bend
    std::vector<uint> neg_splay_ribbon;
    for (int i = 1; i < data->voxels[0]-1; i++) {
        for (int j = 1; j < data->voxels[1]-1; j++) {
            for (int k = 1; k < data->voxels[2]-1; k++) {
                uint id = sub2ind(i, j, k, 0);
                // upper(mag) <= Ssb < lower(mag)
                //if (_widget.knot_interaction_handle.sb_neg_iso_ratio_upper * min_sb <= field_Ssb[id] &&
                //    field_Ssb[id] < min_sb * _widget.knot_interaction_handle.sb_neg_iso_ratio_lower) {
                //    neg_splay_ribbon.emplace_back(id);
                //}

                float local_avg_field_Ssb = field_Ssb[id];
                for (int x : {-1, 1})
                    local_avg_field_Ssb += field_Ssb[sub2ind(i + x, j, k, 0)];
                for (int y : {-1, 1})
                    local_avg_field_Ssb += field_Ssb[sub2ind(i, j + y, k, 0)];
                for (int z : {-1, 1})
                    local_avg_field_Ssb += field_Ssb[sub2ind(i, j, k + z, 0)];

                local_avg_field_Ssb /= 7.;


                if (- _widget.knot_interaction_handle.sb_neg_iso_ratio_upper * stdev < field_Ssb[id]
                    && field_Ssb[id] < - _widget.knot_interaction_handle.sb_neg_iso_ratio_lower * stdev) {
                    neg_splay_ribbon.emplace_back(id);
                }
                valid_fieldf[id] = 0.f;
            }
        }
    }

    // Saturate the ribbon to enhance the contrast
    for (const auto& id : neg_splay_ribbon) {
        // Convert the index to i,j,k space
        int k = id / slice;
        int j = (id - k * slice) / data->voxels[0];
        int i = id - j * data->voxels[0] - slice * k;

        for (int x = -Radius; x <= Radius; x++)
            for (int y = -Radius; y <= Radius; y++)
                for (int z = -Radius; z <= Radius; z++) {
                    // Compute the index to modify
                    uint id_ball = (i + x) + (j + y) * data->voxels[0] + (k + z) * slice;
                    // Make sure defect is in bounds
                    if (id_ball < vol && id_ball >= 0) {
                        //_vortex_line.valid_field[id_ball] = 0;
                        float r2 = x * x + y * y + z * z;
                        valid_fieldf[id_ball] = -1.f;
                    }
                }

    }

    mc.reset_mesh();
    mc.run();

    for (uint i = 0; i < mc.nverts(); ++i) {
        LC::ExtendedMC::Vertex& v = mc.vertices()[i];
        v.x = dr[0] * v.x - cellf[0] * 0.5f;
        v.y = dr[1] * v.y - cellf[1] * 0.5f;
        v.z = dr[2] * v.z - cellf[2] * 0.5f;

        // Normalize normals
        float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
        if (nrm != 0)
        {
            nrm = 1.0 / sqrt(nrm);
            v.nx *= nrm;
            v.ny *= nrm;
            v.nz *= nrm;
        }
    }

    if (mc.nverts() && mc.ntrigs()) {
        unsigned int nVert = mc.nverts();
        unsigned int nInd = mc.ntrigs() * 3;
        LC::Math::IsoVertex* verts = new LC::Math::IsoVertex[nVert];
        unsigned int* indices = new unsigned int[nInd];

        LC::ExtendedMC::Vertex* temp_verts = mc.vertices();
        int* ind_temp = (int*)mc.triangles();

        // Feed data to indices,verts

        for (int i = 0; i < nInd; i++)
            indices[i] = ind_temp[i];

        for (int i = 0; i < nVert; i++) {
            verts[i].position[0] = temp_verts[i].x;
            verts[i].position[1] = temp_verts[i].y;
            verts[i].position[2] = temp_verts[i].z;
            verts[i].normal[0] = temp_verts[i].nx;
            verts[i].normal[1] = temp_verts[i].ny;
            verts[i].normal[2] = temp_verts[i].nz;

        }


        std::vector<MeshLib::PNCVertex<float>> vertices;
        std::vector<MeshLib::Triangle> triangles;
        std::array<float, 4> yellow = { 1.f,1.f,0.f,1.f };

        // Re-interleave data to IsoVertex format
        std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc.nverts()]);
        LC::ExtendedMC::Vertex* iso_v_ptr = mc.vertices();
        uint* iso_ind_ptr = (uint*)mc.triangles();

        for (int vi = 0; vi < mc.nverts(); vi++) {
            iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
        }

        // Smooth the mesh and save it
        //TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(), data->voxels, yellow, false, 0, 0, saveObj, obj_name + std::string("_SB_neg.ply"));
        
        float lambda = 0.33f;
        float mu = -0.34f;
        int smoothIterations = 80;
        TaubinSmoothingGeneral(lambda, mu, smoothIterations, iso_verts.get(), mc.nverts(), iso_ind_ptr,
            mc.ntrigs() * 3, field_nn.get(), data->voxels, yellow, false, 0, 0, saveObj,
            obj_name + std::string("_SB_neg.ply"));


        mc.reset_mesh();
    }

#endif

    // <dog>
    // Compute director preimages
    Eigen::Vector3f cijk[2][2][2];
    Eigen::Vector3f cjk[2][2];
    Eigen::Vector3f ck[2];

    int preim_ct = 0;
    int nPreimages = _widget.knot_interaction_handle.preimages.size();

    // Record the vertices and return
    std::vector<std::tuple<std::vector<LC::Math::IsoVertex>, std::vector<uint>>> heliknoton_list;
    heliknoton_list.reserve(2);

    for (auto preim : _widget.knot_interaction_handle.preimages) {
        float iso = _widget.knot_interaction_handle.preimage_isovalue;
        float th = preim.theta / 180.f * M_PI;
        float ph = preim.phi / 180.f * M_PI;

        Eigen::Vector3f p0(sin(th)* cos(ph), sin(th)* sin(ph), cos(th));

        std::tuple<std::vector<LC::Math::IsoVertex>,std::vector<uint>> p_list;

        auto get_director = [&](int i, int j, int k) {
            // Apply PBCs
            if (i >= data->voxels[0])
                i -= data->voxels[0];
            if (j >= data->voxels[1])
                j -= data->voxels[1];
            if (k >= data->voxels[2])
                k -= data->voxels[2];
            uint32_t v_i = i + j * data->voxels[0] + k * slice;
            return Eigen::Vector3f(field_nn[v_i], field_nn[v_i + vol], field_nn[v_i + 2 * vol]);
        };

        for (int i = 0; i < vox_interp[0]; i++) {
            for (int j = 0; j < vox_interp[1]; j++) {
                for (int k = 0; k < vox_interp[2]; k++) {

                    // interpolant id
                    uint32_t ii = i + j * vox_interp[0] + k * vox_interp[0] * vox_interp[1];

                    // Use trilinear interpolation
                    float x = i * dr_interp[0];
                    float y = j * dr_interp[1];
                    float z = k * dr_interp[2];
                    // Bottom corner points
                    int i_0 = x / dr[0];
                    int j_0 = y / dr[1];
                    int k_0 = z / dr[2];
                    float x_0 = i_0 * dr[0];
                    float y_0 = j_0 * dr[1];
                    float z_0 = k_0 * dr[2];
                    // Interpolation coordinates
                    float x_d = (x - x_0) / dr[0];
                    float y_d = (y - y_0) / dr[1];
                    float z_d = (z - z_0) / dr[2];
                    // Fill the interpolant
                    for (int a = 0; a < 2; a++)
                        for (int b = 0; b < 2; b++)
                            for (int c = 0; c < 2; c++)
                                cijk[a][b][c] = get_director(i_0 + a, j_0 + b, k_0 + c);

                    // Begin interpolation
                    for (int b = 0; b < 2; b++)
                        for (int c = 0; c < 2; c++)
                            cjk[b][c] = (1.f - x_d) * cijk[0][b][c] + x_d * cijk[1][b][c];

                    // Begin interpolation
                    for (int c = 0; c < 2; c++)
                        ck[c] = (1.f - y_d) * cjk[0][c] + y_d * cjk[1][c];

                    Eigen::Vector3f result = (1.f - z_d)* ck[0] + z_d * ck[1];
                    result.normalize();

                    sample_grid[ii] = (result - p0).squaredNorm();
                }
            }
        }

        mc2.reset_mesh();
        mc2.run( iso );

        for (uint i = 0; i < mc2.nverts(); ++i) {
            LC::ExtendedMC::Vertex& v = mc2.vertices()[i];
            v.x = dr_interp[0] * v.x - cellf[0] * 0.5f;
            v.y = dr_interp[1] * v.y - cellf[1] * 0.5f;
            v.z = dr_interp[2] * v.z - cellf[2] * 0.5f;

            // Normalize normals
            float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
            if (nrm != 0)
            {
                nrm = 1.0 / sqrt(nrm);
                v.nx *= nrm;
                v.ny *= nrm;
                v.nz *= nrm;
            }
        }

        if (mc2.nverts() && mc2.ntrigs()) {
            unsigned int nVert = mc2.nverts();
            unsigned int nInd = mc2.ntrigs() * 3;

            std::get<0>(p_list).reserve(nVert);

            LC::ExtendedMC::Vertex* temp_verts = mc2.vertices();
            int* ind_temp = (int*)mc2.triangles();

            std::vector<MeshLib::PNCVertex<float>> vertices;
            std::vector<MeshLib::Triangle> triangles;
            auto col = LC::Imaging::Colors::RungeSphere(th, 2.f * ph);
            std::array<float, 4> col4 = { col[0], col[1], col[2], 1.f };

            // Re-interleave data to IsoVertex format
            std::unique_ptr<LC::Math::IsoVertex[]> iso_verts(new LC::Math::IsoVertex[mc2.nverts()]);
            LC::ExtendedMC::Vertex* iso_v_ptr = mc2.vertices();
            uint* iso_ind_ptr = (uint*)mc2.triangles();

            // Export heliknoton data
            for (int vi = 0; vi < mc2.nverts(); vi++) {
                iso_verts[vi] = extended_to_iso(iso_v_ptr[vi]);
                std::get<0>(p_list).emplace_back(iso_verts[vi]);
            }

            for (int i = 0; i < mc2.ntrigs(); i++) {
                std::get<1>(p_list).emplace_back(iso_ind_ptr[3*i]);
                std::get<1>(p_list).emplace_back(iso_ind_ptr[3*i+1]);
                std::get<1>(p_list).emplace_back(iso_ind_ptr[3*i+2]);
            }

            float lambda = 0.33f;
            float mu = -0.34f;
            int smoothIterations = 80;
            std::string fullName = obj_name + "_preim_" + std::to_string(++preim_ct) + "_of_" + std::to_string(nPreimages) + ".ply";
            TaubinSmoothingGeneral(lambda, mu, smoothIterations, iso_verts.get(), mc2.nverts(), iso_ind_ptr,
                mc2.ntrigs() * 3, field_nn.get(), data->voxels, col4, false, 0, 0, saveObj,
                fullName);
            mc2.reset_mesh();
            // Append to big list
            heliknoton_list.emplace_back(p_list);
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
    Eigen::Vector3f N;

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


#if EXTENDEDMC
    LC::ExtendedMC::MarchingCubes mc;

    mc.clean_all();
    mc.set_resolution(vNew[0], vNew[1], vNew[2]);
    mc.init_all();
#endif

    for (S2Fiber& pimage : _nematicPreimages) {
        // Compute the field (diffmag)
        float theta = pimage.theta / 180.f * M_PI;
        float phi = pimage.phi / 180.f * M_PI;

        N = { sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };

        for (unsigned int i = 0; i < size; i++) {

            if (pimage.pontryagin) {
                field[i] = pow(field_nn[i + 2 * size] - N.z(), 2);
            }
            else
                field[i] = (field_nn[i] - N[0]) * (field_nn[i] - N[0]) +
                (field_nn[i + size] - N[1]) * (field_nn[i + size] - N[1]) +
                (field_nn[i + 2 * size] - N[2]) * (field_nn[i + 2 * size] - N[2]);

        }

        // Count the number of points within the preimage
        unsigned int numPointsFound = 0;
        for (unsigned int i = 0; i < size; i++) if (field[i] < pimage.isoLevel) ++numPointsFound;

        // If points found is greater than half of volume, then invert domain

        bool flip_domain = numPointsFound > size / 2;
        
        // Do the opposite depending on chosen domain inversion
        if (pimage.domain_inversion)
            flip_domain = !flip_domain;

        if (flip_domain) {
            LC_INFO("S2Fiber points found [{0}/{1}] exceeds half of volume: Inverting domain", numPointsFound, size);
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

        // set data

#if EXTENDEDMC
        mc.set_ext_data(field.get());
        mc.run( pimage.isoLevel );
        

        // Rescale positions
        {
            float dx = data->cell_dims[0] / (vNew[0] - 1);
            float dy = data->cell_dims[1] / (vNew[1] - 1);
            float dz = data->cell_dims[2] / (vNew[2] - 1);
            float xmin = -0.5 * data->cell_dims[0];
            float ymin = -0.5 * data->cell_dims[1];
            float zmin = -0.5 * data->cell_dims[2];

            for (uint i = 0; i < mc.nverts(); ++i) {
                LC::ExtendedMC::Vertex& v = mc.vertices()[i];
                v.x = dx * v.x + xmin;
                v.y = dy * v.y + ymin;
                v.z = dz * v.z + zmin;

                // Normalize normals
                float nrm = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
                if (nrm != 0)
                {
                    nrm = 1.0 / sqrt(nrm);
                    v.nx *= nrm;
                    v.ny *= nrm;
                    v.nz *= nrm;
                }
            }
        }

        if (mc.nverts() && mc.ntrigs()) {

#else
        // Generate 

        _isoGenerator.GenerateSurface(field.get(), pimage.isoLevel, vNew, cell, color);

        if (_isoGenerator.isSurfaceValid()) {
#endif

#if EXTENDEDMC
            unsigned int nVert = mc.nverts();
            unsigned int nInd = mc.ntrigs() * 3;
            LC::Math::IsoVertex* verts = new LC::Math::IsoVertex[nVert];
            unsigned int* indices = new unsigned int[nInd];

            LC::ExtendedMC::Vertex* temp_verts = mc.vertices();
            int* ind_temp = (int*)mc.triangles();

            // Feed data to indices,verts

            for (int i = 0; i < nInd; i++)
                indices[i] = ind_temp[i];

            for (int i = 0; i < nVert; i++) {
                verts[i].position[0] = temp_verts[i].x;
                verts[i].position[1] = temp_verts[i].y;
                verts[i].position[2] = temp_verts[i].z;
                verts[i].normal[0] = temp_verts[i].nx;
                verts[i].normal[1] = temp_verts[i].ny;
                verts[i].normal[2] = temp_verts[i].nz;
               
            }

            mc.reset_mesh();

#else
            unsigned int nVert = _isoGenerator.NumSurfaceVertices();
            unsigned int nInd = _isoGenerator.NumSurfaceIndices();
            LC::Math::IsoVertex* verts = _isoGenerator.ReleaseSurfaceVertices();
            unsigned int* indices = _isoGenerator.ReleaseSurfaceIndices();
#endif

            // Subdivide the mesh
            LC::Algorithm::Mesh<float> mesh;
            std::vector<MeshLib::PNCVertex<float>> vertices;
            std::vector<MeshLib::Triangle> triangles;

            mesh.read_pnc_mesh((MeshLib::PNCVertex<float>*)verts, nVert, (MeshLib::Triangle*)indices, nInd / 3);
            LC::Algorithm::loop_subdivision(&mesh, pimage.subdivisions);
            // Export mesh data to vertices and indices (Don't compute normals yet)
            mesh.write_obj(vertices, triangles, false);

            auto smoothing_result = TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(), vNew, color, pimage.normal_inversion, pimage.pontryagin, pimage.subdivisions);

            vertices = std::get<0>(smoothing_result);
            triangles = std::get<1>(smoothing_result);

            //SmoothIsosurface((LC::Math::IsoVertex*)&vertices[0], (unsigned int *)&triangles[0], vertices.size(), triangles.size() * 3, _widget.smoothingIterations, _widget.smoothingValue, _widget.smoothingType);

            LC_INFO("Successfully generated surface (verts = {0}, triangles = {1})", vertices.size(), triangles.size());

            // Fill magnum class with generated surface data
            pimage.surface.Init((LC::Surface::Vertex*)&vertices[0], vertices.size(), &triangles[0][0], triangles.size() * 3, _widget.preimage_translate);

            // Delete data
            delete[] verts;
            delete[] indices;

            // Reset generator
#if !EXTENDEDMC
            _isoGenerator.DeleteSurface();
#endif

            pimage.mesh = pimage.surface.Mesh();

            // Add mesh to the scene
            new LC::Drawable::TransparentNormalDrawable{ *_preimageManipulator, _phongShader, *pimage.mesh, pimage.draw, pimage.cullFaces, _transparentNormalDrawables };

        }
    }

#if EXTENDEDMC
    mc.clean_temps();
#endif

    // Make Vortex Knot isosurface
    std::array<float, 4> color1 = { 1., 1., 0., 1. }, color2 = { 1., 1., 1., 1. };
    std::unique_ptr<float[]> chi_field, shell;
    std::unique_ptr<short[]> valid_field;
    LC::Math::Isosurface<float*, float> gen;
    LC::Math::Isosurface<short*, short> gen2;


    unsigned int ill_defined_chi = 0;
    std::array<float,3> cellf = {(float)data->cell_dims[0],(float)data->cell_dims[1],(float)data->cell_dims[2]};

    if (_widget.chiColorScheme != 2 && _vortexShell) {
        ill_defined_chi = LC::Math::ChiralityField(field_nn.get(), chi_field, vNew, cellf, valid_field, true, _widget.isoLevel);
    }

    if (_widget.generateKnots) {
        // Clean up all knots excluding components
        if (!_processedVortexKnot.empty()) {
            VortexKnot* knot;
            for (auto it = _processedVortexKnot.begin(); it != _processedVortexKnot.end(); it++) {
                knot = &(*it);
                knot->draw = false;
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

#define RBF_METHOD 0

#if RBF_METHOD
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

                LC_INFO("Start point = ({0},{1},{2})", startpoint(0), startpoint(1), startpoint(2));

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

                uint halfCubeSide = max(2. * 3., pow(_vortex_line.knn * 3. / 4. / M_PI, 1. / 3.));
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

                // Compute COM of the loop
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
                    else // Was a duplicate loop
                    {
                        LC_WARN("Duplicate loop ignored.");
                        continue;
                    }

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
#else
        if (_vortex_line.numNodes && _vortex_line.components.size() > 0) {

            // Translate each component into sets of points
            LC::scalar dx = data->cell_dims[0] / (vox[0] - 1);
            LC::scalar dy = data->cell_dims[1] / (vox[1] - 1);
            LC::scalar dz = data->cell_dims[2] / (vox[2] - 1);

            LC::scalar dr = sqrt(dx * dx + dy * dy + dz * dz);
            // Separation distance of two points
            LC::scalar min_point_dist = 2. * _widget.knot_interaction_handle.point_density * dr;


            std::vector<Eigen::Vector3d> distinct_points;
                

            for (const auto& vertex : _vortex_line.vertices) {

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

            // Compute nearest 3 neighbors using knn
            // Convert distinct points into proper data format for knn
            std::unique_ptr<LC::scalar[]> positions(new LC::scalar[distinct_points.size() * 3]);
            std::unique_ptr<std::size_t[]> query_domain(new std::size_t[distinct_points.size()]);
            for (int i = 0; i < distinct_points.size(); i++) {
                // Entire query domain for the knot
                query_domain[i] = i;
                for (int d = 0; d < 3; d++)
                    positions[i + d * distinct_points.size()] = (distinct_points[i])(d);
            }

            int knn = 4;
            LC::Math::Metric<LC::scalar> metric;
            LC::Math::StencilWeights<LC::scalar> derivative = LC::Math::StencilWeightFirstDerivatives<LC::scalar>{};
            metric.Bcs = { false, false, false };
            metric.SetBox(data->cell_dims[0], data->cell_dims[1], data->cell_dims[2]);
            std::unique_ptr<LC::Math::rbf<LC::scalar>> RBF = std::unique_ptr<LC::Math::poly_spline<LC::scalar>>(new LC::Math::poly_spline<LC::scalar>);
            std::unique_ptr<std::size_t[]> neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[distinct_points.size() * knn]);

            // Find nearest neighbors
            LC::Algorithm::knn_c(positions.get(), distinct_points.size(),
                query_domain.get(), distinct_points.size(),
                metric, knn, (LC::scalar*)0, neighbors.get());

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

            LC_INFO("Discarding {0} points", discard_list.size());
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
                    metric, knn, (LC::scalar*)0, neighbors.get());
            }

            // Create a graph with the remaining vertices
            // Add an edge between first neighbor
            Graph graph(parsed_points.size());
            int mapoffset = parsed_points.size();

            for (int i = 0; i < parsed_points.size(); i++) {
                
                // Compute the dist
                std::size_t nb1 = neighbors[1 * mapoffset + i];
                std::size_t nb2 = neighbors[2 * mapoffset + i];
                std::size_t nb3 = neighbors[3 * mapoffset + i];

                Eigen::Vector3d pt = parsed_points[i];
                Eigen::Vector3d pt_nb1 = parsed_points[nb1];
                Eigen::Vector3d pt_nb2 = parsed_points[nb2];
                Eigen::Vector3d pt_nb3 = parsed_points[nb3];

                Eigen::Vector3d t1 = pt - pt_nb1;
                t1.normalize();

                Eigen::Vector3d t2 = pt - pt_nb2;
                t2.normalize();

                Eigen::Vector3d t3 = pt - pt_nb3;
                LC::scalar diff = t3.norm();
                t3.normalize();

                LC::scalar d13 = pow(t1.dot(t3),2);
                LC::scalar d23 = pow(t2.dot(t3),2);


                graph.addEdge(i, nb1);
                graph.addEdge(i, nb2);

                if (diff < 2. * min_point_dist && d13 < 0.5 && d23 < 0.5) {
                    graph.addEdge(i, nb3);
                }
            }
                
            std::list<uint>* adjacencyList = graph.adjacencyList();

            // Export edges and vertices (in lclab2, draw cylinders between each edge connection found)
            
            // Initialize list of visited vertices to zero
            std::map<uint, bool> visited;
            for (int i = 0; i < parsed_points.size(); i++)
                visited.insert({ i,false });





            // Plot the points
            // Vortex knot is valid at this point so generate
            _processedVortexKnot.push_back(VortexKnot{});

            // Export points
            VortexKnot* knot;

            for (auto it = _processedVortexKnot.begin(); it != _processedVortexKnot.end(); it++) {
                knot = &(*it);
            }

            std::function<Vector3(void*, std::size_t)> access = [](void* points, std::size_t i) {
                Vector3 result;
                result[0] = ((Eigen::Vector3d*)points + i)->x();
                result[1] = ((Eigen::Vector3d*)points + i)->y();
                result[2] = ((Eigen::Vector3d*)points + i)->z();
                return result;
            };

            {
                std::list<VortexKnot>::iterator iter = _vortexKnot.begin();
                knot->knotColor = iter->knotColor;
            }

            knot->points.polyRadius = 0.5 * min_point_dist;
            knot->points.Init(&distinct_points[0].x(), access, distinct_points.size());

            // TODO
            // knot->cylinders....

            // Identify components
            auto components = graph.connectedComponents(3);

            std::vector<std::array<Eigen::Vector3d, 2>> edges;


            //for (int i = 0; i < parsed_points.size(); i++) {
            //    for (auto j : adjacencyList[i]) {
            //        if (!visited[j]) {
            //            edges.push_back({ parsed_points[i],parsed_points[j] });
            //        }
            //    }
            //    visited[i] = true;
            //}
            std::vector<int> nedges_per_component(components.size());
            for (int i = 0; i < components.size(); i++) {
                nedges_per_component[i] = 0;
                for (const auto& ci : components[i]) {
                    for (auto j : adjacencyList[ci]) {
                        if (!visited[j]) {
                            edges.push_back({ parsed_points[ci],parsed_points[j] });
                            ++nedges_per_component[i];
                        }
                    }
                    visited[ci] = true;
                }
            }

            // Export parsed points
            //std::ofstream ofile("D:/dev/lclab2/data/knot/test/points.bin", std::ios::out | std::ios::binary);

            //if (ofile.is_open()) {
            
            //    int numComponents = components.size();
            //    ofile.write((char*)&numComponents, sizeof(int));
                // Each component consists of the number of edges, followed by the edgelist
            //    int offset = 0;
            //    for (int i = 0; i < numComponents; i++) {
            //        ofile.write((char*)&nedges_per_component[i], sizeof(int));
            //        ofile.write((char*)&edges[offset], sizeof(Eigen::Vector3d) * 2 * nedges_per_component[i]);
                    // Increment the offset
            //        offset += nedges_per_component[i];
            //    }
            //    ofile.close();
            //}


        }
#endif
        
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
            LC_INFO("S2Fiber points found [{0}/{1}] exceeds half of volume: Inverting domain", helicalPointsFound, vol);
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

            //SmoothIsosurface(verts, indices, nVert, nInd, 1, 200.f, 2);
            auto smoothing_result = TaubinSmoothing(verts, nVert, indices, nInd, field_nn.get(),
                vNew, {1.f,1.f,1.f,1.f}, false, false, 0);

            auto vertices = std::get<0>(smoothing_result);
            auto triangles = std::get<1>(smoothing_result);

            // Go through vertices and color according to vortex shell tilt
            // Note needs to include preimage_translate since embedded volumes can be chosen!

            // Reduced quaternion space SU2/Q8
            Eigen::Matrix3d I3, A1, A2, A3;
            I3 = Eigen::Matrix3d::Identity();
            A1 = Eigen::DiagonalMatrix<LC::scalar, 3, 3>(1., -1., -1.);
            A2 = Eigen::DiagonalMatrix<LC::scalar, 3, 3>(-1., 1., -1.);
            A3 = Eigen::DiagonalMatrix<LC::scalar, 3, 3>(-1., -1., 1.);

            // Find the min and max chi_x
            float chi_x_min = 1.;
            float chi_x_max = -1.;
            if (_widget.chiColorScheme != 2) {
                for (int i = 0; i < nVert; i++) {
                    int xx = (vertices[i].position[0] / data->cell_dims[0] + 0.5) * (vNew[0] - 1);
                    int yy = (vertices[i].position[1] / data->cell_dims[1] + 0.5) * (vNew[1] - 1);
                    int zz = (vertices[i].position[2] / data->cell_dims[2] + 0.5) * (vNew[2] - 1);
                    unsigned int idx = xx + yy * vNew[0] + zz * slc;
                    float chi_x = chi_field[idx];
                    chi_x_min = chi_x_min > chi_x ? chi_x : chi_x_min;
                    chi_x_max = chi_x_max < chi_x ? chi_x : chi_x_max;
                }
            }

            float dx = data->cell_dims[0] / (vNew[0] - 1);
            float dy = data->cell_dims[1] / (vNew[1] - 1);
            float dz = data->cell_dims[2] / (vNew[2] - 1);

            for (int i = 0; i < nVert; i++) {
                // extract indices from position
                int xx = (vertices[i].position[0] / data->cell_dims[0] + 0.5) * (vNew[0] - 1);
                int yy = (vertices[i].position[1] / data->cell_dims[1] + 0.5) * (vNew[1] - 1);
                int zz = (vertices[i].position[2] / data->cell_dims[2] + 0.5) * (vNew[2] - 1);


                // Get chi
                Color3 tilt_color;

                // Chosen color scheme
                // 0 == tilt (theta)
                // 1 == tilt (chi_x)
                // 2 == quaternion coloring
                if (_widget.chiColorScheme != 2) {

                    // (xx,yy,zz) is the bottom vertex of the voxel containing this (vx,vy,vz) position
                    // Use the 8 corners to linearly interpolate what chi_x, chi_y, and chi_z is at the center
                    
                    // Center the vertex into cube of [0,1]^3
                    float xd = (vertices[i].position[0] - xx * dx) / dx;
                    float yd = (vertices[i].position[1] - yy * dy) / dy;
                    float zd = (vertices[i].position[2] - zz * dz) / dz;


                    // Data mapping
                    // 0 -> (0,0,0)
                    // 1 -> (1,0,0)
                    // 2 -> (0,1,0)
                    // 3 -> (1,1,0)
                    // 4 -> (0,0,1)
                    // 5 -> (1,0,1)
                    // 6 -> (0,1,1)
                    // 7 -> (1,1,1)

                    std::array<float, 8> c_ijk_x;
                    std::array<float, 8> c_ijk_y;
                    std::array<float, 8> c_ijk_z;

                    std::array<float, 4> c_jk_x;
                    std::array<float, 4> c_jk_y;
                    std::array<float, 4> c_jk_z;

                    std::array<float, 2> c_k_x;
                    std::array<float, 2> c_k_y;
                    std::array<float, 2> c_k_z;

                    for (int id_z = 0; id_z < 1; id_z++)
                        for (int id_y = 0; id_y < 1; id_y++)
                            for (int id_x = 0; id_x < 1; id_x++)
                            {
                                int id = id_x + 2 * id_y + 4 * id_z;
                                unsigned int idx = (xx + id_x) + (yy + id_y) * vNew[0] + (zz + id_z) * slc;

                                if (idx >= vol) {
                                    c_ijk_x[id] = 0.f;
                                    c_ijk_y[id] = 0.f;
                                    c_ijk_z[id] = 0.f;
                                }
                                else {
                                    c_ijk_x[id] = chi_field[idx];
                                    c_ijk_y[id] = chi_field[idx + size];
                                    c_ijk_z[id] = chi_field[idx + 2 * size];
                                }
                            }


                    // Begin interpolation
                    for (int id_z = 0; id_z < 1; id_z++)
                        for (int id_y = 0; id_y < 1; id_y++)
                        {
                            c_jk_x[id_y + 2 * id_z] = c_ijk_x[2 * id_y + 4 * id_z] * (1.f - xd) + c_ijk_x[1 + 2 * id_y + 4 * id_z] * xd;
                            c_jk_y[id_y + 2 * id_z] = c_ijk_y[2 * id_y + 4 * id_z] * (1.f - xd) + c_ijk_y[1 + 2 * id_y + 4 * id_z] * xd;
                            c_jk_z[id_y + 2 * id_z] = c_ijk_z[2 * id_y + 4 * id_z] * (1.f - xd) + c_ijk_z[1 + 2 * id_y + 4 * id_z] * xd;
                        }

                    for (int id_z = 0; id_z < 1; id_z++) {
                        c_k_x[id_z] = c_jk_x[2 * id_z] * (1.f - yd) + c_jk_x[1 + 2 * id_z] * yd;
                        c_k_y[id_z] = c_jk_y[2 * id_z] * (1.f - yd) + c_jk_y[1 + 2 * id_z] * yd;
                        c_k_z[id_z] = c_jk_z[2 * id_z] * (1.f - yd) + c_jk_z[1 + 2 * id_z] * yd;
                    }


                    float chi_x = c_k_x[0] * (1.f - zd) + c_k_x[1] * zd;
                    float chi_y = c_k_y[0] * (1.f - zd) + c_k_y[1] * zd;
                    float chi_z = c_k_z[0] * (1.f - zd) + c_k_z[1] * zd;

                    // Normalize
                    LC::Math::Normalize(chi_x, chi_y, chi_z);

                    // Flip chi -> -chi if chi_z does not align with far-field (chosen to be +zhat)
                    if (chi_z < 0) {
                        chi_x = -chi_x;
                        chi_y = -chi_y;
                        chi_z = -chi_z;
                    }


                    // Get angle in xy plane
                    if (_widget.chiColorScheme == 0)
                    {
                        float theta = acos(chi_z);
                        float phi = 2 * atan2(chi_y, chi_x);
                        if (phi < 0.0f) phi = phi + 2 * M_PI;

                        tilt_color = LC::Imaging::Colors::RungeSphere(M_PI/2., phi);
                        //tilt_color = pow(cos(angle2), 2) * purple + pow(sin(angle2),2) * tilt_color;
                    }
                    else if (_widget.chiColorScheme == 1)
                    {
                        Color3 neutral = Color3(.29f, 0.f, 0.51f); // Indigo
                        Color3 cyan = Color3::cyan();
                        Color3 yellow = Color3::yellow();

                        // Map chi_x into range [0,1]
                        float zbar;
                        if (chi_x_max != chi_x_min)
                            zbar = (chi_x - chi_x_min) / (chi_x_max - chi_x_min);
                        else
                            zbar = 1.0;

                        Color3 result;

                        // Map negative values to yellow: zbar in [0, 0.5]
                        if (chi_x < 0.) result = (1. - zbar / .5) * yellow + zbar / .5 * neutral;
                        // Map positive values to cyan: zbar in [0.5, 1]
                        else result = (1. - (zbar - 0.5) / .5) * neutral + (zbar - 0.5) / .5 * cyan;

                        tilt_color = result;
                    }
                }
                else {
                    unsigned int idx = xx + yy * data->voxels[0] + zz * data->voxels[0] * data->voxels[1];
                    tilt_color = Color3::blue();
                    Eigen::Quaterniond q{ 0.,0.,0.,1. };
                    if (_quaternion_field)
                        q = _quaternion_field[idx];

#define SU2modQ8 0
#if SU2modQ8
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
#else
                    // Color mapping
                    // 0 -> negative qx : green
                    // 1 -> positive qx : pink
                    // 2 -> negative qy : cyan
                    // 3 -> positive qy : yellow
                    // 4 -> negative qz : dark blue
                    // 5 -> positive qz : red
                    int col_sel = q.x() > 0. ? 1 : 0;
                    float qmax = abs(q.x());

                    if (qmax < abs(q.y()))
                    {
                        qmax = abs(q.y());
                        col_sel = q.y() > 0. ? 3 : 2;
                    }
                    
                    if (qmax < abs(q.z())) {
                        qmax = abs(q.z());
                        col_sel = q.z() > 0. ? 5 : 4;
                    }

                    switch (col_sel) {
                        case 0:
                            tilt_color = Color3::green();
                            break;
                        case 1:
                            tilt_color = Color3(1.f, .41f, .71f); // pink
                            break;
                        case 2:
                            tilt_color = Color3::cyan();
                            break;
                        case 3:
                            tilt_color = Color3::yellow();
                            break;
                        case 4:
                            tilt_color = Color3::blue();
                            break;
                        case 5:
                            tilt_color = Color3::red();
                            break;
                        default:
                            LC_WARN("Quaternion color cases failed");
                    }

#endif
                }

                // Apply the color
                for (int d = 0; d < 3; d++) {

                    vertices[i].color[d] = tilt_color[d];

                    if (_vortexShell->invertNormals)
                        vertices[i].normal[d] *= -1;
                }

            }

            LC_INFO("Successfully generated surface (verts = {0}, indices = {1})", nVert, nInd);

            // Fill magnum class with generated surface data
            _vortexShell->surface.Init((LC::Surface::Vertex*)&vertices[0], vertices.size(), &triangles[0][0], triangles.size() * 3, _widget.preimage_translate);

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
            LC_INFO("S2Fiber points found [{0}/{1}] exceeds half of volume: Inverting domain", domainPoints, baryVol);
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
            LC_INFO("S2Fiber points found [{0}/{1}] exceeds half of volume: Inverting domain", domainPoints, baryVol);
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

std::tuple<std::vector<Eigen::Vector3d>,std::vector<int>> exportVortexKnot(const std::vector<MeshLib::PNCVertex<float>>& vertices, const std::vector<MeshLib::Triangle>& triangles, LC::scalar min_point_dist, const std::string &fname) {

    std::vector<Eigen::Vector3d> points;
    std::vector<int> npoints_per_component;
    std::map<uint, Face> unvisited_map;
    // Split triangles in components
    uint nTriangles = triangles.size();

    // Note that verts currently refers to verts and not global indices
    // Therefore, convert indices data to global index space through vertex positions
    for (int tri = 0; tri < nTriangles; tri++) {
        uint i1 = triangles[tri][0];
        uint i2 = triangles[tri][1];
        uint i3 = triangles[tri][2];

        // Store in map
        unvisited_map.insert({ tri, Face(i1, i2, i3) });
    }

    // Find components through the vertices
    auto vortex_components = find_all_components_graph(unvisited_map, vertices.size(), 50);

    LC_INFO("Vortex components detected - {0}", vortex_components.size());

    if (vortex_components.empty())
        return {};

    std::ofstream ofile(fname.c_str(), std::ios::out | std::ios::binary);

    // Convert each component into a set of points
    std::vector<std::vector<Eigen::Vector3d>> distinct_component_points;

    for (int c = 0; c < vortex_components.size(); c++) {
        std::vector<Eigen::Vector3d> pointset;
        for (auto id : vortex_components[c]) {
            auto vertex = vertices[id];
            Eigen::Vector3d potential_pos(vertex.position[0], vertex.position[1], vertex.position[2]);
            // Check if position is too close to existing points
            bool valid_pt = true;
            for (const auto& pt : pointset) {
                LC::scalar dist = (potential_pos - pt).norm();
                if (dist < min_point_dist) {
                    valid_pt = false;
                    break;
                }
            }

            // Add to the distinct points
            if (valid_pt) {
                pointset.emplace_back(potential_pos);
            }
        }

        LC_INFO("Component {1}: {0} distinct points detected", pointset.size(),c);
        
        if (pointset.size() < 5) {
            LC_WARN("Not enough points detected, discarding component");
            continue;
        }
        distinct_component_points.push_back(pointset);
    }

    int knn = 5;
    LC::Math::Metric<LC::scalar> metric;
    metric.Bcs = { false, false, false };
    std::unique_ptr<LC::Math::rbf<LC::scalar>> RBF = std::unique_ptr<LC::Math::poly_spline<LC::scalar>>(new LC::Math::poly_spline<LC::scalar>);

    // For each component...
    for (auto& distinct_points : distinct_component_points) {
        // See how many distinct points
        // Compute nearest 2 neighbors using knn
        // Convert distinct points into proper data format for knn
        std::unique_ptr<LC::scalar[]> positions(new LC::scalar[distinct_points.size() * 3]);
        std::unique_ptr<std::size_t[]> query_domain(new std::size_t[distinct_points.size()]);
        for (int i = 0; i < distinct_points.size(); i++) {
            // Entire query domain for the knot
            query_domain[i] = i;
            for (int d = 0; d < 3; d++)
                positions[i + d * distinct_points.size()] = (distinct_points[i])(d);
        }

        std::unique_ptr<std::size_t[]> neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[distinct_points.size() * knn]);

        // Find nearest neighbors
        LC::Algorithm::knn_c(positions.get(), distinct_points.size(),
            query_domain.get(), distinct_points.size(),
            metric, knn, (LC::scalar*)0, neighbors.get());

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
                metric, knn, (LC::scalar*)0, neighbors.get());
        }

        // Create a graph with the remaining vertices
        // Add an edge between first neighbor

        Graph graph(parsed_points.size());
        int mapoffset = parsed_points.size();

        for (int i = 0; i < parsed_points.size(); i++) {

            std::vector<std::pair<uint,Eigen::Vector3d>> tangents;

            for (int k = 1; k < 3; k++) {
                std::size_t nb_k = neighbors[k * mapoffset + i];

                Eigen::Vector3d pt = parsed_points[i];
                Eigen::Vector3d pt_nbk = parsed_points[nb_k];

                Eigen::Vector3d tk = pt - pt_nbk;
                LC::scalar diff = tk.norm();
                tk.normalize();

                if (diff < 2. * min_point_dist) {
                    tangents.push_back({ nb_k, tk });
                }
            }

            // if there are only two tangents, insert both
            if (tangents.size() == 2) {
                graph.addEdge(i, tangents[0].first);
                graph.addEdge(i, tangents[1].first);
            }
            else if (tangents.size() > 2) {
                // Determine which connections are the best

                // Triple intersection point
                if (tangents.size() == 3) {
                    // Cannot determine yet so add all of the points
                    graph.addEdge(i, tangents[0].first);
                    graph.addEdge(i, tangents[1].first);
                    graph.addEdge(i, tangents[2].first);
                }

                // Two knots intersecting each other
                if (tangents.size() == 4) {
                    // Choose the two connections with the dot product of largest magnitude
                    // There are 4 choose 2 connections -> 6 total possibilities
                    std::array<LC::scalar, 6> prods;

                    /*
                    * Index: connection
                     0: 0--1
                     1: 1--2
                     2: 2--3
                     3: 1--3
                     4: 0--2
                     5: 0--3
                    */
                    prods[0] = abs(tangents[0].second.dot(tangents[1].second));
                    prods[1] = abs(tangents[1].second.dot(tangents[2].second));
                    prods[2] = abs(tangents[2].second.dot(tangents[3].second));
                    prods[3] = abs(tangents[1].second.dot(tangents[3].second));
                    prods[4] = abs(tangents[0].second.dot(tangents[2].second));
                    prods[5] = abs(tangents[0].second.dot(tangents[3].second));

                    int max_id = 0;

                    for (int d = 1; d < 6; d++) {
                        if (prods[d] > prods[max_id])
                            max_id = d;
                    }

                    int c1, c2, c3, c4;

                    if (max_id == 0) {
                        c1 = 0;
                        c2 = 1;
                        c3 = 2;
                        c4 = 3;
                    }
                    else if (max_id == 1) {
                        c1 = 1;
                        c2 = 2;
                        c3 = 3;
                        c4 = 0;
                    }
                    else if (max_id == 2) {
                        c1 = 2;
                        c2 = 3;
                        c3 = 1;
                        c4 = 0;
                    }
                    else if (max_id == 3) {
                        c1 = 1;
                        c2 = 3;
                        c3 = 2;
                        c4 = 0;
                    }
                    else if (max_id == 4) {
                        c1 = 0;
                        c2 = 2;
                        c3 = 1;
                        c4 = 3;
                    }
                    else if (max_id == 5) {
                        c1 = 0;
                        c2 = 3;
                        c3 = 2;
                        c4 = 1;
                    }
                    graph.addEdge(tangents[c1].first, tangents[c2].first);
                    graph.addEdge(tangents[c3].first, tangents[c4].first);
                }

            }

        }

        std::list<uint>* adjacencyList = graph.adjacencyList();

        // Export edges and vertices (in lclab2, draw cylinders between each edge connection found)

        // Initialize list of visited vertices to zero
        std::map<uint, bool> visited;
        for (int i = 0; i < parsed_points.size(); i++)
            visited.insert({ i,false });

        // Automatically sort indices
        auto components = graph.connectedComponents(3);

        for (int i = 0; i < components.size(); i++) {
            npoints_per_component.push_back(0);
            int end = npoints_per_component.size() - 1;
            npoints_per_component[end] = 0;
            for (const auto& ci : components[i]) {

                points.emplace_back(parsed_points[ci]);
                ++npoints_per_component[end];
            }
        }
    
    }

    if (points.size()) {
        int numComponents = npoints_per_component.size();
        if (ofile.is_open()) {

            ofile.write((char*)&numComponents, sizeof(int));
            // Each component consists of the number of edges, followed by the edgelist
            int offset = 0;
            for (int i = 0; i < numComponents; i++) {
                ofile.write((char*)&npoints_per_component[i], sizeof(int));
                ofile.write((char*)&points[offset], sizeof(Eigen::Vector3d) * npoints_per_component[i]);
                // Increment the offset
                offset += npoints_per_component[i];
            }
            ofile.close();
        }

            // Save secondary format for knotplot as txt format
        std::string base_fname = fname.substr(0, fname.size() - 4);

        std::ofstream ofile_knotplot((base_fname+".txt").c_str(), std::ios::out);
        if (ofile_knotplot.is_open()) {
            int offset = 0;
            for (const auto& ni : npoints_per_component) {
                for (int i = 0; i < ni; i++) {
                    ofile_knotplot << (float)points[offset + i].x() * 20.f << " "
                        << (float)points[offset + i].y() * 20.f << " "
                        << (float)points[offset + i].z() * 20.f << std::endl;
                }
                offset += ni;
                ofile_knotplot << std::endl;
            }
        }

    }
    return { points, npoints_per_component };
}


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