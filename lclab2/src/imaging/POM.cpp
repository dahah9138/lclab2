#include "POM.h"

namespace LC { namespace Imaging { namespace UniformGrid {

    // Analytical CIE 1931 XYZ fit (Wyman et al. 2013) -> Linear sRGB weights for wavelength in nm
    static std::array<scalar, 3> spectralToLinearRGB(scalar lambda_nm) {
        auto asymGauss = [](scalar x, scalar mu, scalar t1, scalar t2) -> scalar {
            scalar t = (x - mu) * ((x < mu) ? t1 : t2);
            return std::exp(-0.5 * t * t);
            };

        scalar X = 1.056 * asymGauss(lambda_nm, 599.8, 0.0264, 0.0323)
            + 0.362 * asymGauss(lambda_nm, 442.0, 0.0624, 0.0374)
            - 0.065 * asymGauss(lambda_nm, 501.1, 0.0490, 0.0382);

        scalar Y = 0.821 * asymGauss(lambda_nm, 568.8, 0.0213, 0.0247)
            + 0.286 * asymGauss(lambda_nm, 530.9, 0.0613, 0.0322);

        scalar Z = 1.217 * asymGauss(lambda_nm, 437.0, 0.0845, 0.0278)
            + 0.681 * asymGauss(lambda_nm, 459.0, 0.0385, 0.0725);

        scalar r = 3.2406 * X - 1.5372 * Y - 0.4986 * Z;
        scalar g = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
        scalar b = 0.0557 * X - 0.2040 * Y + 1.0570 * Z;

        return { std::max<scalar>(0.0, r), std::max<scalar>(0.0, g), std::max<scalar>(0.0, b) };
    }

    void POM::Compute(scalar* nn, const std::array<int, 3>& voxels, void* CData, ColorDataFunc colorFunc, const float& alpha) {
        if (lights.empty()) return;

        const std::complex<scalar> ii(0, 1);
        const scalar polAngleRad = polarizerAngle * M_PI / 180.0;
        const scalar anaAngleRad = polAngleRad + M_PI / 2.0;
        const std::size_t numLights = lights.size();

        std::array<std::size_t, 3> rawSlices = {
            static_cast<std::size_t>(voxels[0]),
            static_cast<std::size_t>(voxels[0]) * voxels[1],
            static_cast<std::size_t>(voxels[0]) * voxels[1] * voxels[2]
        };
        std::array<std::size_t, 3> cumVoxProdSlice = {
            rawSlices[indexOrder[0]],
            rawSlices[indexOrder[1]],
            rawSlices[indexOrder[2]]
        };

        std::vector<scalar> k0_dz(numLights);
        std::vector<std::complex<scalar>> eid0_vec(numLights);
        std::vector<std::array<scalar, 3>> rgbWeights(numLights);
        std::array<scalar, 3> weightNorm = { 0.0, 0.0, 0.0 };

        for (std::size_t l = 0; l < numLights; ++l) {
            const scalar lambda_m = std::max<scalar>(1e-9, lights[l].wavelength * 1e-9);
            k0_dz[l] = (2.0 * M_PI / lambda_m) * dz;
            eid0_vec[l] = std::exp(ii * (k0_dz[l] * n0));

            auto w = spectralToLinearRGB(lights[l].wavelength);
            for (int c = 0; c < 3; ++c) {
                rgbWeights[l][c] = w[c] * lights[l].intensity;
                weightNorm[c] += w[c];
            }
        }

        // Normalize spectral matching sum so pure white at I=1 stays near 1.0
        const scalar maxChannelNorm = std::max({ weightNorm[0], weightNorm[1], weightNorm[2], scalar(1e-6) });
        for (std::size_t l = 0; l < numLights; ++l) {
            for (int c = 0; c < 3; ++c) {
                rgbWeights[l][c] /= maxChannelNorm;
            }
        }

        const int scan_layer_depth = std::clamp(static_cast<int>(z_scan_ratio * voxels[2]), 0, voxels[2]);

        using PolarizationState = std::vector<Eigen::Vector2cd>;
        std::vector<PolarizationState> Eo_states;

        if (polarizer) {
            Eo_states.push_back(PolarizationState(numLights, { std::cos(polAngleRad), std::sin(polAngleRad) }));
        }
        else {
            Eo_states.push_back(PolarizationState(numLights, { 1.0, 0.0 }));
            Eo_states.push_back(PolarizationState(numLights, { 0.0, 1.0 }));
        }

        const scalar invNumStates = 1.0 / static_cast<scalar>(Eo_states.size());

        if (additional_layers > 0 && dop > 0.0f) {
            const scalar nodes_in_layer = voxels[2] / (dop * 2.0);
            const scalar nodes_in_layers = nodes_in_layer * additional_layers;
            const scalar dPhi = (nodes_in_layer > 1.0) ? (1.0 / (nodes_in_layer - 1.0)) : 0.0;
            Eigen::Matrix2cd M;

            for (std::size_t l = 0; l < numLights; ++l) {
                const std::complex<scalar> eidE = std::exp(ii * (k0_dz[l] * ne));
                const std::complex<scalar> eid0 = eid0_vec[l];

                for (int t = 0; t < static_cast<int>(nodes_in_layers) - 1; ++t) {
                    const scalar phi = M_PI * dPhi * t + M_PI / 2.0;
                    const scalar cp = std::cos(phi);
                    const scalar sp = std::sin(phi);

                    M(0, 0) = cp * cp * eidE + sp * sp * eid0;
                    M(0, 1) = sp * cp * (eidE - eid0);
                    M(1, 0) = M(0, 1);
                    M(1, 1) = sp * sp * eidE + cp * cp * eid0;

                    for (auto& Eo_plane : Eo_states) {
                        Eo_plane[l] = M * Eo_plane[l];
                    }
                }
            }
        }

        // Precompute waveplate matrices outside the (i, j) pixel loop
        std::vector<Eigen::Matrix2cd> wpMatrices(numLights, Eigen::Matrix2cd::Identity());
        if (waveplate != Waveplate::None) {
            const scalar full_wave_530 = 530e-9 / (1.55338 - 1.54425);
            const scalar thick = (waveplate == Waveplate::Full530nm) ? full_wave_530 : (0.25 * full_wave_530);

            for (std::size_t l = 0; l < numLights; ++l) {
                const scalar lambda_m = std::max<scalar>(1e-9, lights[l].wavelength * 1e-9);
                const scalar de_wp = 2.0 * M_PI / lambda_m * thick * 1.55338;
                const scalar do_ep = 2.0 * M_PI / lambda_m * thick * 1.54425;

                const std::complex<scalar> expE = std::exp(ii * de_wp);
                const std::complex<scalar> expO = std::exp(ii * do_ep);

                wpMatrices[l](0, 0) = 0.5 * (expE + expO);
                wpMatrices[l](0, 1) = 0.5 * (expE - expO);
                wpMatrices[l](1, 0) = wpMatrices[l](0, 1);
                wpMatrices[l](1, 1) = wpMatrices[l](0, 0);
            }
        }

        const scalar ca = std::cos(anaAngleRad);
        const scalar sa = std::sin(anaAngleRad);

        for (int i = 0; i < voxels[0]; ++i) {
            for (int j = 0; j < voxels[1]; ++j) {
                auto Eo_vector = Eo_states;
                Eigen::Matrix2cd M;

                for (int k = 0; k < scan_layer_depth; ++k) {
                    scalar nx = nn[dir2ind(i, j, k, 0, cumVoxProdSlice)];
                    scalar ny = nn[dir2ind(i, j, k, 1, cumVoxProdSlice)];
                    scalar nz = nn[dir2ind(i, j, k, 2, cumVoxProdSlice)];

                    const scalar n_norm = std::sqrt(nx * nx + ny * ny + nz * nz);
                    if (n_norm > 0.0) {
                        nx /= n_norm; ny /= n_norm; nz /= n_norm;
                    }

                    const scalar ct = nz;
                    const scalar st2 = nx * nx + ny * ny;
                    const scalar st = std::sqrt(st2);
                    const scalar cp = (st > 1e-12) ? (nx / st) : 1.0;
                    const scalar sp = (st > 1e-12) ? (ny / st) : 0.0;
                    const scalar ne_th = ne * n0 / std::sqrt(n0 * n0 * st2 + ne * ne * ct * ct);

                    for (std::size_t l = 0; l < numLights; ++l) {
                        const std::complex<scalar> eidE = std::exp(ii * (k0_dz[l] * ne_th));
                        const std::complex<scalar> eid0 = eid0_vec[l];

                        M(0, 0) = cp * cp * eidE + sp * sp * eid0;
                        M(0, 1) = sp * cp * (eidE - eid0);
                        M(1, 0) = M(0, 1);
                        M(1, 1) = sp * sp * eidE + cp * cp * eid0;

                        for (auto& Eo : Eo_vector) {
                            Eo[l] = M * Eo[l];
                        }
                    }
                }

                if (waveplate != Waveplate::None) {
                    for (std::size_t l = 0; l < numLights; ++l) {
                        for (auto& Eo : Eo_vector) {
                            Eo[l] = wpMatrices[l] * Eo[l];
                        }
                    }
                }

                std::array<scalar, 3> accumRGB = { 0.0, 0.0, 0.0 };

                for (std::size_t l = 0; l < numLights; ++l) {
                    scalar I_lambda = 0.0;
                    for (const auto& Eo : Eo_vector) {
                        if (analyzer) {
                            std::complex<scalar> E_out = Eo[l](0) * ca + Eo[l](1) * sa;
                            I_lambda += std::norm(E_out);
                        }
                        else {
                            I_lambda += Eo[l].squaredNorm();
                        }
                    }
                    I_lambda *= invNumStates;

                    accumRGB[0] += rgbWeights[l][0] * I_lambda;
                    accumRGB[1] += rgbWeights[l][1] * I_lambda;
                    accumRGB[2] += rgbWeights[l][2] * I_lambda;
                }

                std::array<float, 4> color = {
                    static_cast<float>(std::pow(std::max<scalar>(0.0, accumRGB[0]), gamma)),
                    static_cast<float>(std::pow(std::max<scalar>(0.0, accumRGB[1]), gamma)),
                    static_cast<float>(std::pow(std::max<scalar>(0.0, accumRGB[2]), gamma)),
                    alpha
                };

                colorFunc(CData, color, cross2ind(i, j, cumVoxProdSlice[0]));
            }
        }
    }

}}}