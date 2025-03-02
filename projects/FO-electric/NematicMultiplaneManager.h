#include <lclab2.h>


struct NematicMultiplaneManager {
	using FOFDSolver = LC::FrankOseen::Electric::FOFDSolver;
	using Dataset = FOFDSolver::Dataset;
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