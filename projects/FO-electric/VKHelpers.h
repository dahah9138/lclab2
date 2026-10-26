#ifndef VORTEXHELPER_H
#define VORTEXHELPER_H

#include "find_components.h"
#include <math.h>

namespace VKH {
    struct GridMetrics {
        std::array<int, 3> dims;
        std::array<float, 3> cell;
        std::array<float, 3> dr;
        uint slice;
        uint vol;

        GridMetrics(const std::array<int, 3>& vox, const std::array<float, 3>& cell_dims, int mult = 1) {
            for (int d = 0; d < 3; ++d) {
                dims[d] = vox[d] * mult;
                cell[d] = cell_dims[d];
                dr[d] = cell[d] / static_cast<float>(dims[d] - 1);
            }
            slice = dims[0] * dims[1];
            vol = slice * dims[2];
        }

        inline uint to1D(int i, int j, int k) const {
            return i + dims[0] * j + slice * k;
        }

        inline void to3D(uint id, int& i, int& j, int& k) const {
            k = id / slice;
            uint rem = id % slice;
            j = rem / dims[0];
            i = rem % dims[0];
        }

        inline bool inBounds(int i, int j, int k) const {
            return i >= 0 && i < dims[0] &&
                j >= 0 && j < dims[1] &&
                k >= 0 && k < dims[2];
        }

        inline bool inInterior(int i, int j, int k, int margin = 2) const {
            return i >= margin && i < dims[0] - margin &&
                j >= margin && j < dims[1] - margin &&
                k >= margin && k < dims[2] - margin;
        }
    };

    // Forward declare the helper
    float evalSplayBendAtNode(const GridMetrics& grid, const float* field_Q, int i, int j, int k);

    struct SmoothedMesh {
        std::vector<MeshLib::PNCVertex<float>> vertices;
        std::vector<MeshLib::Triangle> triangles;

        bool empty() const { return vertices.empty() || triangles.empty(); }
    };

    LC::Math::IsoVertex toIsoVertex(const LC::ExtendedMC::Vertex& v) {
        LC::Math::IsoVertex viso;
        viso.position[0] = v.x; viso.position[1] = v.y; viso.position[2] = v.z;
        viso.normal[0] = v.nx;  viso.normal[1] = v.ny;  viso.normal[2] = v.nz;
        return viso;
    }

    void rescaleAndNormalizeMC(LC::ExtendedMC::MarchingCubes& mc, const GridMetrics& grid) {
        for (uint i = 0; i < mc.nverts(); ++i) {
            LC::ExtendedMC::Vertex& v = mc.vertices()[i];
            v.x = grid.dr[0] * v.x - 0.5f * grid.cell[0];
            v.y = grid.dr[1] * v.y - 0.5f * grid.cell[1];
            v.z = grid.dr[2] * v.z - 0.5f * grid.cell[2];

            float nrm2 = v.nx * v.nx + v.ny * v.ny + v.nz * v.nz;
            if (nrm2 > 0.f) {
                float inv_nrm = 1.0f / std::sqrt(nrm2);
                v.nx *= inv_nrm;
                v.ny *= inv_nrm;
                v.nz *= inv_nrm;
            }
        }
    }

    void rescaleAndNormalizeMesh(SmoothedMesh& mesh, const GridMetrics& grid) {
        for (auto& v : mesh.vertices) {
            for (int d = 0; d < 3; ++d) {
                v.position[d] = grid.dr[d] * v.position[d] - 0.5f * grid.cell[d];
            }
            float nrm2 = v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] + v.normal[2] * v.normal[2];
            if (nrm2 > 0.f) {
                float inv_nrm = 1.0f / std::sqrt(nrm2);
                v.normal[0] *= inv_nrm;
                v.normal[1] *= inv_nrm;
                v.normal[2] *= inv_nrm;
            }
        }
    }

    // Update signature and tau normalization inside the existing helper:
    void computeOrderParameterAndFilter(
        const GridMetrics& grid,
        const float* field_nn,
        const float* chi_field,
        float iso_value,
        short* valid_field,
        float* valid_fieldf,
        float* field_S,
        bool normalize_tau = false)
    {
        for (uint ii = 0; ii < grid.vol; ++ii) {
            valid_fieldf[ii] = 0.f;
            field_S[ii] = 0.f;

            if (!valid_field[ii]) continue;

            int i, j, k;
            grid.to3D(ii, i, j, k);

            Eigen::Vector3d chi0(chi_field[ii], chi_field[ii + grid.vol], chi_field[ii + 2 * grid.vol]);
            Eigen::Vector3d nn0(field_nn[ii], field_nn[ii + grid.vol], field_nn[ii + 2 * grid.vol]);
            Eigen::Vector3d tau0 = chi0.cross(nn0);
            if (normalize_tau && tau0.squaredNorm() > 0.) {
                tau0.normalize();
            }

            LC::scalar avgcos2_chi = 0.;
            LC::scalar avgcos2_tau = 0.;
            int points_in_bds = 0;

            for (int z = -1; z <= 1; ++z) {
                for (int y = -1; y <= 1; ++y) {
                    for (int x = -1; x <= 1; ++x) {
                        if (x == 0 && y == 0 && z == 0) continue;
                        if (!grid.inBounds(i + x, j + y, k + z)) continue;

                        uint ip = grid.to1D(i + x, j + y, k + z);
                        points_in_bds++;

                        Eigen::Vector3d chi_nb(chi_field[ip], chi_field[ip + grid.vol], chi_field[ip + 2 * grid.vol]);
                        Eigen::Vector3d nn_nb(field_nn[ip], field_nn[ip + grid.vol], field_nn[ip + 2 * grid.vol]);
                        Eigen::Vector3d tau_nb = chi_nb.cross(nn_nb);
                        if (normalize_tau && tau_nb.squaredNorm() > 0.) {
                            tau_nb.normalize();
                        }

                        double dot_chi = chi0.dot(chi_nb);
                        double dot_tau = tau0.dot(tau_nb);
                        avgcos2_chi += dot_chi * dot_chi;
                        avgcos2_tau += dot_tau * dot_tau;
                    }
                }
            }

            if (points_in_bds > 0) {
                LC::scalar Sop_chi = 0.5 * (3. * avgcos2_chi / points_in_bds - 1.);
                LC::scalar Sop_tau = 0.5 * (3. * avgcos2_tau / points_in_bds - 1.);
                field_S[ii] = static_cast<float>(0.5 * (Sop_chi + Sop_tau));
            }
        }

        for (uint ii = 0; ii < grid.vol; ++ii) {
            if (field_S[ii] < iso_value) {
                valid_fieldf[ii] = -1.f;
                valid_field[ii] = 0;
            }
        }
    }


    float computeBaseGridSplayBend(
        const GridMetrics& base_grid,
        const float* valid_fieldf,
        const float* field_Q,
        int sb_rad,
        float* field_Ssb)
    {
        std::vector<uint8_t> fat_mask(base_grid.vol, 0);

        for (uint id = 0; id < base_grid.vol; ++id) {
            if (valid_fieldf[id] >= 0.f) continue;

            int i, j, k;
            base_grid.to3D(id, i, j, k);

            for (int z = -sb_rad; z <= sb_rad; ++z) {
                for (int y = -sb_rad; y <= sb_rad; ++y) {
                    for (int x = -sb_rad; x <= sb_rad; ++x) {
                        if (base_grid.inBounds(i + x, j + y, k + z)) {
                            fat_mask[base_grid.to1D(i + x, j + y, k + z)] = 1;
                        }
                    }
                }
            }
        }

        std::fill_n(field_Ssb, base_grid.vol, 0.f);
        double sum = 0.0;

        for (int k = 2; k < base_grid.dims[2] - 2; ++k) {
            for (int j = 2; j < base_grid.dims[1] - 2; ++j) {
                for (int i = 2; i < base_grid.dims[0] - 2; ++i) {
                    uint id = base_grid.to1D(i, j, k);
                    if (fat_mask[id]) {
                        float val = evalSplayBendAtNode(base_grid, field_Q, i, j, k);
                        field_Ssb[id] = val;
                        sum += val;
                    }
                }
            }
        }

        float mean = static_cast<float>(sum / static_cast<double>(base_grid.vol));
        double var_sum = 0.0;
        for (uint i = 0; i < base_grid.vol; ++i) {
            double diff = field_Ssb[i] - mean;
            var_sum += diff * diff;
        }

        return (base_grid.vol > 0) ? static_cast<float>(std::sqrt(var_sum / base_grid.vol)) : 0.f;
    }

    void populatePreimageSampleGrid(
        const GridMetrics& base_grid,
        const GridMetrics& interp_grid,
        const float* field_nn,
        const Eigen::Vector3f& p0,
        float* sample_grid)
    {
        auto get_director = [&](int i, int j, int k) {
            if (i >= base_grid.dims[0]) i -= base_grid.dims[0];
            if (j >= base_grid.dims[1]) j -= base_grid.dims[1];
            if (k >= base_grid.dims[2]) k -= base_grid.dims[2];
            uint32_t v_i = base_grid.to1D(i, j, k);
            return Eigen::Vector3f(
                field_nn[v_i],
                field_nn[v_i + base_grid.vol],
                field_nn[v_i + 2 * base_grid.vol]
            );
            };

        for (int k = 0; k < interp_grid.dims[2]; ++k) {
            float z = k * interp_grid.dr[2];
            int k0 = static_cast<int>(z / base_grid.dr[2]);
            float zd = (z - k0 * base_grid.dr[2]) / base_grid.dr[2];

            for (int j = 0; j < interp_grid.dims[1]; ++j) {
                float y = j * interp_grid.dr[1];
                int j0 = static_cast<int>(y / base_grid.dr[1]);
                float yd = (y - j0 * base_grid.dr[1]) / base_grid.dr[1];

                for (int i = 0; i < interp_grid.dims[0]; ++i) {
                    float x = i * interp_grid.dr[0];
                    int i0 = static_cast<int>(x / base_grid.dr[0]);
                    float xd = (x - i0 * base_grid.dr[0]) / base_grid.dr[0];

                    Eigen::Vector3f cijk[2][2][2];
                    for (int a = 0; a < 2; ++a)
                        for (int b = 0; b < 2; ++b)
                            for (int c = 0; c < 2; ++c)
                                cijk[a][b][c] = get_director(i0 + a, j0 + b, k0 + c);

                    Eigen::Vector3f cjk[2][2];
                    for (int b = 0; b < 2; ++b)
                        for (int c = 0; c < 2; ++c)
                            cjk[b][c] = (1.f - xd) * cijk[0][b][c] + xd * cijk[1][b][c];

                    Eigen::Vector3f ck[2];
                    for (int c = 0; c < 2; ++c)
                        ck[c] = (1.f - yd) * cjk[0][c] + yd * cjk[1][c];

                    Eigen::Vector3f result = (1.f - zd) * ck[0] + zd * ck[1];
                    result.normalize();

                    sample_grid[interp_grid.to1D(i, j, k)] = (result - p0).squaredNorm();
                }
            }
        }
    }

    std::vector<std::vector<uint>> extractConnectedComponents(const SmoothedMesh& mesh) {
        std::map<uint, Face> unvisited_faces;
        uint nTriangles = static_cast<uint>(mesh.triangles.size());

        for (uint tri = 0; tri < nTriangles; ++tri) {
            unvisited_faces.insert({ tri, Face(mesh.triangles[tri][0], mesh.triangles[tri][1], mesh.triangles[tri][2]) });
        }

        try {
            return find_all_components_graph(unvisited_faces, static_cast<uint>(mesh.vertices.size()), 15);
        }
        catch (const std::runtime_error& e) {
            LC_INFO("Error: {0}", e.what());
            return {};
        }
    }

    void stampComponentToGrid(
        const std::vector<uint>& component,
        const std::vector<MeshLib::PNCVertex<float>>& vertices,
        const GridMetrics& grid_interp,
        int radius,
        float* sample_grid)
    {
        std::fill_n(sample_grid, grid_interp.vol, 0.1f);

        for (uint ci : component) {
            const auto& pos = vertices[ci].position;
            int i = static_cast<int>((pos[0] + 0.5f * grid_interp.cell[0]) / grid_interp.cell[0] * (grid_interp.dims[0] - 1));
            int j = static_cast<int>((pos[1] + 0.5f * grid_interp.cell[1]) / grid_interp.cell[1] * (grid_interp.dims[1] - 1));
            int k = static_cast<int>((pos[2] + 0.5f * grid_interp.cell[2]) / grid_interp.cell[2] * (grid_interp.dims[2] - 1));

            for (int z = -radius; z <= radius; ++z) {
                for (int y = -radius; y <= radius; ++y) {
                    for (int x = -radius; x <= radius; ++x) {
                        if (!grid_interp.inBounds(i + x, j + y, k + z)) continue;

                        uint id_ball = grid_interp.to1D(i + x, j + y, k + z);
                        float r2 = static_cast<float>(x * x + y * y + z * z);
                        sample_grid[id_ball] += (r2 == 0.f) ? -1.f : (-1.f / r2);
                    }
                }
            }
        }
    }

    void colorVerticesByDirector(
        std::vector<MeshLib::PNCVertex<float>>& vertices,
        const GridMetrics& base_grid,
        const GridMetrics& interp_grid,
        const float* field_nn,
        std::vector<uint8_t>& upsampled_mask,
        std::vector<uint>& upsampled_indices)
    {
        for (auto& v : vertices) {
            int x = static_cast<int>((v.position[0] / base_grid.cell[0] + 0.5f) * (base_grid.dims[0] - 1));
            int y = static_cast<int>((v.position[1] / base_grid.cell[1] + 0.5f) * (base_grid.dims[1] - 1));
            int z = static_cast<int>((v.position[2] / base_grid.cell[2] + 0.5f) * (base_grid.dims[2] - 1));

            int xi = static_cast<int>((v.position[0] / interp_grid.cell[0] + 0.5f) * (interp_grid.dims[0] - 1));
            int yi = static_cast<int>((v.position[1] / interp_grid.cell[1] + 0.5f) * (interp_grid.dims[1] - 1));
            int zi = static_cast<int>((v.position[2] / interp_grid.cell[2] + 0.5f) * (interp_grid.dims[2] - 1));

            if (interp_grid.inBounds(xi, yi, zi)) {
                uint idi = interp_grid.to1D(xi, yi, zi);
                if (!upsampled_mask[idi]) {
                    upsampled_mask[idi] = 1;
                    upsampled_indices.push_back(idi);
                }
            }

            if (base_grid.inBounds(x, y, z)) {
                uint id = base_grid.to1D(x, y, z);
                float theta = std::acos(field_nn[id + 2 * base_grid.vol]);
                float phi = std::atan2(field_nn[id + base_grid.vol], field_nn[id]);
                Magnum::Color3 dcol = LC::Imaging::Colors::RungeSphere(theta, phi);
                v.color[0] = dcol[0];
                v.color[1] = dcol[1];
                v.color[2] = dcol[2];
                v.color[3] = 1.f;
            }
        }
    }

    void colorVerticesByDirector(
        std::vector<MeshLib::PNCVertex<float>>& vertices,
        const GridMetrics& base_grid,
        const float* field_nn)
    {
        for (auto& v : vertices) {
            int x = static_cast<int>((v.position[0] / base_grid.cell[0] + 0.5f) * (base_grid.dims[0] - 1));
            int y = static_cast<int>((v.position[1] / base_grid.cell[1] + 0.5f) * (base_grid.dims[1] - 1));
            int z = static_cast<int>((v.position[2] / base_grid.cell[2] + 0.5f) * (base_grid.dims[2] - 1));

            if (base_grid.inBounds(x, y, z)) {
                uint id = base_grid.to1D(x, y, z);
                float theta = std::acos(field_nn[id + 2 * base_grid.vol]);
                float phi = std::atan2(field_nn[id + base_grid.vol], field_nn[id]);
                Magnum::Color3 dcol = LC::Imaging::Colors::RungeSphere(theta, phi);
                v.color[0] = dcol[0];
                v.color[1] = dcol[1];
                v.color[2] = dcol[2];
                v.color[3] = 1.f;
            }
        }
    }

    void computeChiQTensor(const GridMetrics& grid, const float* chi_field, float* field_Q) {
        for (uint ii = 0; ii < grid.vol; ++ii) {
            double cx = chi_field[ii];
            double cy = chi_field[ii + grid.vol];
            double cz = chi_field[ii + 2 * grid.vol];

            field_Q[ii] = static_cast<float>(0.5 * (3.0 * cx * cx - 1.0));
            field_Q[ii + grid.vol] = static_cast<float>(1.5 * cx * cy);
            field_Q[ii + 2 * grid.vol] = static_cast<float>(1.5 * cx * cz);
            field_Q[ii + 3 * grid.vol] = static_cast<float>(0.5 * (3.0 * cy * cy - 1.0));
            field_Q[ii + 4 * grid.vol] = static_cast<float>(1.5 * cy * cz);
        }
    }

    float evalSplayBendAtNode(const GridMetrics& grid, const float* field_Q, int i, int j, int k) {
        auto q = [&](int x, int y, int z, int comp) {
            return field_Q[grid.to1D(x, y, z) + comp * grid.vol];
            };

        float qa002 = (q(i, j, k + 1, 0) + q(i, j, k - 1, 0) - 2.f * q(i, j, k, 0)) / (grid.dr[2] * grid.dr[2]);
        float qa200 = (q(i + 1, j, k, 0) + q(i - 1, j, k, 0) - 2.f * q(i, j, k, 0)) / (grid.dr[0] * grid.dr[0]);
        float qd002 = (q(i, j, k + 1, 3) + q(i, j, k - 1, 3) - 2.f * q(i, j, k, 3)) / (grid.dr[2] * grid.dr[2]);
        float qd020 = (q(i, j + 1, k, 3) + q(i, j - 1, k, 3) - 2.f * q(i, j, k, 3)) / (grid.dr[1] * grid.dr[1]);

        float qb110 = (q(i + 1, j + 1, k, 1) + q(i - 1, j - 1, k, 1) - q(i - 1, j + 1, k, 1) - q(i + 1, j - 1, k, 1)) / (4.f * grid.dr[0] * grid.dr[1]);
        float qc101 = (q(i + 1, j, k + 1, 2) + q(i - 1, j, k - 1, 2) - q(i - 1, j, k + 1, 2) - q(i + 1, j, k - 1, 2)) / (4.f * grid.dr[0] * grid.dr[2]);
        float qe011 = (q(i, j + 1, k + 1, 4) + q(i, j - 1, k - 1, 4) - q(i, j - 1, k + 1, 4) - q(i, j + 1, k - 1, 4)) / (4.f * grid.dr[1] * grid.dr[2]);

        return -qa002 + qa200 + 2.f * (qb110 + qc101) - qd002 + qd020 + 2.f * qe011;
    }

    std::vector<uint> buildFatVortexRegion(
        const GridMetrics& interp_grid,
        const std::vector<uint>& upsampled_indices,
        const std::vector<uint8_t>& upsampled_mask,
        int sb_rad)
    {
        std::vector<uint8_t> fat_mask(interp_grid.vol, 0);
        std::vector<uint> fat_list;

        for (uint id : upsampled_indices) {
            int i, j, k;
            interp_grid.to3D(id, i, j, k);

            for (int z = -sb_rad; z <= sb_rad; ++z) {
                for (int y = -sb_rad; y <= sb_rad; ++y) {
                    for (int x = -sb_rad; x <= sb_rad; ++x) {
                        if (!interp_grid.inBounds(i + x, j + y, k + z)) continue;

                        uint id_ball = interp_grid.to1D(i + x, j + y, k + z);
                        if (!upsampled_mask[id_ball] && !fat_mask[id_ball]) {
                            fat_mask[id_ball] = 1;
                            fat_list.push_back(id_ball);
                        }
                    }
                }
            }
        }
        return fat_list;
    }

    struct SplayBendStats {
        float max_sb = 0.f;
        float min_sb = 0.f;
        float tc_pos = 0.f;
        float tc_neg = 0.f;
    };

    SplayBendStats computeInterpolatedSplayBend(
        const GridMetrics& base_grid,
        const GridMetrics& interp_grid,
        const float* field_Q,
        const std::vector<uint>& fat_vortex_list,
        float* field_Ssb)
    {
        std::fill_n(field_Ssb, interp_grid.vol, 0.f);

        for (uint id : fat_vortex_list) {
            int i, j, k;
            interp_grid.to3D(id, i, j, k);

            if (!interp_grid.inInterior(i, j, k, 2)) continue;

            float x = i * interp_grid.dr[0];
            float y = j * interp_grid.dr[1];
            float z = k * interp_grid.dr[2];

            int i0 = static_cast<int>(x / base_grid.dr[0]);
            int j0 = static_cast<int>(y / base_grid.dr[1]);
            int k0 = static_cast<int>(z / base_grid.dr[2]);

            if (!base_grid.inInterior(i0, j0, k0, 1) || !base_grid.inInterior(i0 + 1, j0 + 1, k0 + 1, 1)) continue;

            float xd = (x - i0 * base_grid.dr[0]) / base_grid.dr[0];
            float yd = (y - j0 * base_grid.dr[1]) / base_grid.dr[1];
            float zd = (z - k0 * base_grid.dr[2]) / base_grid.dr[2];

            float cijk[2][2][2];
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    for (int c = 0; c < 2; ++c)
                        cijk[a][b][c] = evalSplayBendAtNode(base_grid, field_Q, i0 + a, j0 + b, k0 + c);

            float cjk[2][2];
            for (int b = 0; b < 2; ++b)
                for (int c = 0; c < 2; ++c)
                    cjk[b][c] = (1.f - xd) * cijk[0][b][c] + xd * cijk[1][b][c];

            float ck[2];
            for (int c = 0; c < 2; ++c)
                ck[c] = (1.f - yd) * cjk[0][c] + yd * cjk[1][c];

            field_Ssb[id] = (1.f - zd) * ck[0] + zd * ck[1];
        }

        SplayBendStats stats;
        float mean_pos = 0.f, mean_neg = 0.f;
        int n_pos = 0, n_neg = 0;

        for (uint id = 0; id < interp_grid.vol; ++id) {
            float val = field_Ssb[id];
            if (val > 0.f) {
                mean_pos += val;
                n_pos++;
                if (val > stats.max_sb) stats.max_sb = val;
            }
            else if (val < 0.f) {
                mean_neg += val;
                n_neg++;
                if (val < stats.min_sb) stats.min_sb = val;
            }
        }

        if (n_pos > 0) mean_pos /= static_cast<float>(n_pos);
        if (n_neg > 0) mean_neg /= static_cast<float>(n_neg);

        constexpr float alpha = 0.1f;
        stats.tc_pos = (stats.max_sb != 0.f) ? (mean_pos / stats.max_sb) * alpha : 0.f;
        stats.tc_neg = (stats.min_sb != 0.f) ? (mean_neg / stats.min_sb) * alpha : 0.f;

        return stats;
    }

    void populateSplayBendSampleGrid(
        const GridMetrics& interp_grid,
        const std::vector<uint>& fat_vortex_list,
        const float* field_Ssb,
        float threshold_ratio,
        float extremum_sb,
        float* sample_grid)
    {
        std::fill_n(sample_grid, interp_grid.vol, 0.f);
        if (extremum_sb == 0.f) return;

        for (uint id : fat_vortex_list) {
            int i, j, k;
            interp_grid.to3D(id, i, j, k);
            if (!interp_grid.inInterior(i, j, k, 1)) continue;

            float local_avg = field_Ssb[id]
                + field_Ssb[interp_grid.to1D(i - 1, j, k)] + field_Ssb[interp_grid.to1D(i + 1, j, k)]
                + field_Ssb[interp_grid.to1D(i, j - 1, k)] + field_Ssb[interp_grid.to1D(i, j + 1, k)]
                + field_Ssb[interp_grid.to1D(i, j, k - 1)] + field_Ssb[interp_grid.to1D(i, j, k + 1)];
            local_avg /= 7.f;

            if (threshold_ratio <= local_avg / extremum_sb) {
                sample_grid[id] = -1.f;
            }
        }
    }

    inline void applyDomainInversion(float* field, uint vol, float iso_level, bool invert_flag = false) {
        uint count = 0;
        for (uint i = 0; i < vol; ++i) {
            if (field[i] < iso_level) ++count;
        }

        bool flip = (count > vol / 2);
        if (invert_flag) flip = !flip;

        if (flip) {
            LC_INFO("S2Fiber points found [{0}/{1}] exceeds half of volume: Inverting domain", count, vol);
            for (uint i = 0; i < vol; ++i) {
                field[i] = (field[i] > iso_level) ? 0.0f : 10.0f;
            }
        }
    }

    inline void buildProcessedKnotGraph(
        const std::vector<MeshLib::PNCVertex<float>>& raw_vertices,
        LC::scalar min_point_dist,
        const std::array<float, 3>& cellf,
        std::vector<Eigen::Vector3d>& distinct_points)
    {
        for (const auto& vertex : raw_vertices) {
            Eigen::Vector3d pos(vertex.position[0], vertex.position[1], vertex.position[2]);
            bool valid_pt = true;
            for (const auto& pt : distinct_points) {
                if ((pos - pt).norm() < min_point_dist) {
                    valid_pt = false;
                    break;
                }
            }
            if (valid_pt) {
                distinct_points.emplace_back(pos);
            }
        }

        LC_INFO("{0} distinct points detected", distinct_points.size());
        if (distinct_points.size() < 4) return;

        size_t n_pts = distinct_points.size();
        auto positions = std::make_unique<LC::scalar[]>(n_pts * 3);
        auto query_domain = std::make_unique<std::size_t[]>(n_pts);

        for (size_t i = 0; i < n_pts; ++i) {
            query_domain[i] = i;
            for (int d = 0; d < 3; ++d) {
                positions[i + d * n_pts] = distinct_points[i](d);
            }
        }

        constexpr int knn = 4;
        LC::Math::Metric<LC::scalar> metric;
        metric.Bcs = { false, false, false };
        metric.SetBox(cellf[0], cellf[1], cellf[2]);
        auto neighbors = std::make_unique<std::size_t[]>(n_pts * knn);

        LC::Algorithm::knn_c(
            positions.get(), n_pts, query_domain.get(), n_pts,
            metric, knn, static_cast<LC::scalar*>(nullptr), neighbors.get()
        );

        std::vector<Eigen::Vector3d> parsed_points;
        parsed_points.reserve(n_pts);
        size_t discarded = 0;

        for (size_t i = 0; i < n_pts; ++i) {
            std::size_t n1 = neighbors[n_pts + i];
            if ((distinct_points[i] - distinct_points[n1]).norm() <= 2.0 * min_point_dist) {
                parsed_points.push_back(distinct_points[i]);
            }
            else {
                ++discarded;
            }
        }

        LC_INFO("Discarding {0} points", discarded);
        if (parsed_points.size() < 4) return;

        if (discarded > 0) {
            size_t n_parsed = parsed_points.size();
            for (size_t i = 0; i < n_parsed; ++i) {
                for (int d = 0; d < 3; ++d) {
                    positions[i + d * n_parsed] = parsed_points[i](d);
                }
            }
            LC::Algorithm::knn_c(
                positions.get(), n_parsed, query_domain.get(), n_parsed,
                metric, knn, static_cast<LC::scalar*>(nullptr), neighbors.get()
            );
        }

        size_t n_parsed = parsed_points.size();
        Graph graph(n_parsed);

        for (size_t i = 0; i < n_parsed; ++i) {
            std::size_t nb1 = neighbors[1 * n_parsed + i];
            std::size_t nb2 = neighbors[2 * n_parsed + i];
            std::size_t nb3 = neighbors[3 * n_parsed + i];

            Eigen::Vector3d t1 = (parsed_points[i] - parsed_points[nb1]).normalized();
            Eigen::Vector3d t2 = (parsed_points[i] - parsed_points[nb2]).normalized();
            Eigen::Vector3d t3 = parsed_points[i] - parsed_points[nb3];
            LC::scalar diff = t3.norm();
            if (diff > 0.) t3 /= diff;

            LC::scalar d13 = t1.dot(t3);
            LC::scalar d23 = t2.dot(t3);

            graph.addEdge(i, nb1);
            graph.addEdge(i, nb2);

            if (diff < 2.0 * min_point_dist && (d13 * d13) < 0.5 && (d23 * d23) < 0.5) {
                graph.addEdge(i, nb3);
            }
        }

        auto components = graph.connectedComponents(3);
        std::list<uint>* adjacencyList = graph.adjacencyList();
        std::vector<uint8_t> visited(n_parsed, 0);
        std::vector<std::array<Eigen::Vector3d, 2>> edges;
        std::vector<int> nedges_per_component(components.size(), 0);

        for (size_t i = 0; i < components.size(); ++i) {
            for (uint ci : components[i]) {
                for (uint j : adjacencyList[ci]) {
                    if (!visited[j]) {
                        edges.push_back({ parsed_points[ci], parsed_points[j] });
                        ++nedges_per_component[i];
                    }
                }
                visited[ci] = 1;
            }
        }
    }

    inline void computeHelicalShellField(
        const GridMetrics& sub_grid,
        float* field_nn,
        float query_value,
        float iso_level,
        float* shell)
    {
        Eigen::Matrix3d handedness0 = Eigen::Matrix3d::Zero();
        handedness0(2, 2) = -1.0;

        std::array<int, 3> vNew = sub_grid.dims;
        std::array<float, 3> cellf = sub_grid.cell;

        for (int k = 0; k < sub_grid.dims[2]; ++k) {
            for (int j = 0; j < sub_grid.dims[1]; ++j) {
                for (int i = 0; i < sub_grid.dims[0]; ++i) {
                    uint idx = sub_grid.to1D(i, j, k);
                    Eigen::Matrix3d handedness = Eigen::Matrix3d::Zero();

                    if (sub_grid.inInterior(i, j, k, 2)) {
                        handedness = (0.5 / M_PI) * LC::Math::HandednessTensor(i, j, k, field_nn, vNew, cellf);
                    }
                    else {
                        handedness(2, 2) = -1.0;
                    }

                    double R2 = 0.5 * (handedness - handedness0).squaredNorm();
                    shell[idx] = static_cast<float>(std::abs(R2 - query_value));
                }
            }
        }

        applyDomainInversion(shell, sub_grid.vol, iso_level, false);
    }

    inline void colorVortexShellVertices(
        std::vector<MeshLib::PNCVertex<float>>& vertices,
        const GridMetrics& sub_grid,
        const GridMetrics& base_grid,
        const float* chi_field,
        const Eigen::Quaterniond* quaternion_field,
        int color_scheme,
        bool invert_normals)
    {
        float chi_x_min = 1.f;
        float chi_x_max = -1.f;

        if (color_scheme != 2 && chi_field) {
            for (const auto& v : vertices) {
                int xx = static_cast<int>((v.position[0] / sub_grid.cell[0] + 0.5f) * (sub_grid.dims[0] - 1));
                int yy = static_cast<int>((v.position[1] / sub_grid.cell[1] + 0.5f) * (sub_grid.dims[1] - 1));
                int zz = static_cast<int>((v.position[2] / sub_grid.cell[2] + 0.5f) * (sub_grid.dims[2] - 1));
                if (sub_grid.inBounds(xx, yy, zz)) {
                    float cx = chi_field[sub_grid.to1D(xx, yy, zz)];
                    chi_x_min = (std::min)(chi_x_min, cx);
                    chi_x_max = (std::max)(chi_x_max, cx);
                }
            }
        }

        for (auto& v : vertices) {
            float fx = (v.position[0] / sub_grid.cell[0] + 0.5f) * (sub_grid.dims[0] - 1);
            float fy = (v.position[1] / sub_grid.cell[1] + 0.5f) * (sub_grid.dims[1] - 1);
            float fz = (v.position[2] / sub_grid.cell[2] + 0.5f) * (sub_grid.dims[2] - 1);

            int xx = static_cast<int>(fx);
            int yy = static_cast<int>(fy);
            int zz = static_cast<int>(fz);

            Color3 tilt_color = Color3::blue();

            if (color_scheme != 2 && chi_field) {
                float xd = fx - xx;
                float yd = fy - yy;
                float zd = fz - zz;

                Eigen::Vector3f c_ijk[2][2][2];
                for (int dz = 0; dz < 2; ++dz) {
                    for (int dy = 0; dy < 2; ++dy) {
                        for (int dx = 0; dx < 2; ++dx) {
                            if (sub_grid.inBounds(xx + dx, yy + dy, zz + dz)) {
                                uint idx = sub_grid.to1D(xx + dx, yy + dy, zz + dz);
                                c_ijk[dx][dy][dz] = Eigen::Vector3f(
                                    chi_field[idx],
                                    chi_field[idx + sub_grid.vol],
                                    chi_field[idx + 2 * sub_grid.vol]
                                );
                            }
                            else {
                                c_ijk[dx][dy][dz] = Eigen::Vector3f::Zero();
                            }
                        }
                    }
                }

                Eigen::Vector3f c_jk[2][2];
                for (int dy = 0; dy < 2; ++dy)
                    for (int dz = 0; dz < 2; ++dz)
                        c_jk[dy][dz] = (1.f - xd) * c_ijk[0][dy][dz] + xd * c_ijk[1][dy][dz];

                Eigen::Vector3f c_k[2];
                for (int dz = 0; dz < 2; ++dz)
                    c_k[dz] = (1.f - yd) * c_jk[0][dz] + yd * c_jk[1][dz];

                Eigen::Vector3f chi = (1.f - zd) * c_k[0] + zd * c_k[1];
                float chi_x = chi.x(), chi_y = chi.y(), chi_z = chi.z();
                LC::Math::Normalize(chi_x, chi_y, chi_z);

                if (chi_z < 0.f) {
                    chi_x = -chi_x;
                    chi_y = -chi_y;
                    chi_z = -chi_z;
                }

                if (color_scheme == 0) {
                    float phi = 2.f * std::atan2(chi_y, chi_x);
                    if (phi < 0.f) phi += 2.f * static_cast<float>(M_PI);
                    tilt_color = LC::Imaging::Colors::RungeSphere(static_cast<float>(M_PI) / 2.f, phi);
                }
                else if (color_scheme == 1) {
                    Color3 neutral(0.29f, 0.f, 0.51f);
                    Color3 cyan = Color3::cyan();
                    Color3 yellow = Color3::yellow();

                    float zbar = (chi_x_max != chi_x_min) ? (chi_x - chi_x_min) / (chi_x_max - chi_x_min) : 1.0f;
                    if (chi_x < 0.f) {
                        tilt_color = (1.f - zbar / 0.5f) * yellow + (zbar / 0.5f) * neutral;
                    }
                    else {
                        tilt_color = (1.f - (zbar - 0.5f) / 0.5f) * neutral + ((zbar - 0.5f) / 0.5f) * cyan;
                    }
                }
            }
            else {
                Eigen::Quaterniond q(0., 0., 0., 1.);
                if (quaternion_field && base_grid.inBounds(xx, yy, zz)) {
                    q = quaternion_field[base_grid.to1D(xx, yy, zz)];
                }

                int col_sel = q.x() > 0. ? 1 : 0;
                double qmax = std::abs(q.x());

                if (qmax < std::abs(q.y())) {
                    qmax = std::abs(q.y());
                    col_sel = q.y() > 0. ? 3 : 2;
                }
                if (qmax < std::abs(q.z())) {
                    col_sel = q.z() > 0. ? 5 : 4;
                }

                switch (col_sel) {
                case 0: tilt_color = Color3::green(); break;
                case 1: tilt_color = Color3(1.f, 0.41f, 0.71f); break;
                case 2: tilt_color = Color3::cyan(); break;
                case 3: tilt_color = Color3::yellow(); break;
                case 4: tilt_color = Color3::blue(); break;
                case 5: tilt_color = Color3::red(); break;
                default: LC_WARN("Quaternion color cases failed"); break;
                }
            }

            for (int d = 0; d < 3; ++d) {
                v.color[d] = tilt_color[d];
                if (invert_normals) v.normal[d] *= -1.f;
            }
        }
    }

    inline void colorVerticesByMaxPion(
        LC::Math::IsoVertex* verts,
        uint nVert,
        const GridMetrics& grid,
        const Eigen::Quaterniond* quaternion_field)
    {
        if (!quaternion_field) return;

        for (uint i = 0; i < nVert; ++i) {
            int xx = static_cast<int>((verts[i].position[0] / grid.cell[0] + 0.5f) * (grid.dims[0] - 1));
            int yy = static_cast<int>((verts[i].position[1] / grid.cell[1] + 0.5f) * (grid.dims[1] - 1));
            int zz = static_cast<int>((verts[i].position[2] / grid.cell[2] + 0.5f) * (grid.dims[2] - 1));

            if (!grid.inBounds(xx, yy, zz)) continue;

            const auto& q = quaternion_field[grid.to1D(xx, yy, zz)];
            LC::scalar max_comp = std::abs(q.x());
            Color3 tilt_color = Color3::blue();

            if (max_comp < std::abs(q.y())) {
                tilt_color = Color3::red();
                max_comp = std::abs(q.y());
            }
            if (max_comp < std::abs(q.z())) {
                tilt_color = Color3::green();
            }

            for (int d = 0; d < 3; ++d) {
                verts[i].color[d] = tilt_color[d];
            }
        }
    }

}

#endif