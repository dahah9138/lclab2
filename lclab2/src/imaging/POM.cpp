#include "POM.h"

namespace LC { namespace Imaging { namespace UniformGrid {

void POM::Compute(scalar *nn, const std::array<int, 3> &voxels, void *CData, ColorDataFunc colorFunc, const float &alpha) {

    // Jones matrices
    Eigen::Matrix2cd M, m;

    scalar th = polarizers * polarizerAngle * M_PI / 180.0;
    const std::complex<scalar> ii(0, 1);


    std::array<std::size_t, 3> cumVoxProdSlice = { voxels[0], voxels[1], voxels[2] };
    std::array<double, 3> lambda = {lightRGB[0] * 1e-9, lightRGB[1] * 1e-9, lightRGB[2] * 1e-9 };

    // slices
    for (int i = 1; i < 3; i++)
        cumVoxProdSlice[i] *= cumVoxProdSlice[i - 1];

    // Permute with indexOrder
    for (int i = 0; i < 3; i++)
        cumVoxProdSlice[i] = cumVoxProdSlice[indexOrder[i]];



    scalar theta, phi, nx, ny, nz;
        for (int i = 0; i < voxels[0]; i++) {
            for (int j = 0; j < voxels[1]; j++) {

                // Polarization states for rgb colors
                Eigen::Vector2cd Eo[] = { {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0} };

                for (int k = 0; k < voxels[2]; k++) {

                    nx = nn[dir2ind(i, j, k, 0, cumVoxProdSlice)];
                    ny = nn[dir2ind(i, j, k, 1, cumVoxProdSlice)];
                    nz = nn[dir2ind(i, j, k, 2, cumVoxProdSlice)];

                    // Compute theta and phi
                    theta = M_PI / 2.0 - atan2(nz, sqrt(nx * nx + ny * ny));
                    phi = atan2(ny, nx);

                    const scalar ct = cos(theta);
                    const scalar st = sin(theta);
                    const scalar cp = cos(phi);
                    const scalar sp = sin(phi);


                    for (int rgb = 0; rgb < 3; rgb++) {
                        const scalar delta0 = 2.0 * M_PI / lambda[rgb] * dz * n0;
                        const scalar ne_th = ne * n0 / sqrt(pow(n0 * st, 2.0) + pow(ne * ct, 2.0));
                        const scalar deltaE = 2.0 * M_PI / lambda[rgb] * dz * ne_th;

                        const std::complex<scalar> eidE = std::exp(ii * deltaE);
                        const std::complex<scalar> eid0 = std::exp(ii * delta0);

                        M(0, 0) = cp * cp * eidE + sp * sp * eid0;
                        M(0, 1) = sp * cp * (eidE - eid0);
                        M(1, 0) = M(0, 1);
                        M(1, 1) = sp * sp * eidE + cp * cp * eid0;

                        Eo[rgb] = M * Eo[rgb];
                    }


                }

                // Additional waveplate
                if (waveplate == Waveplate::Full530nm) {

                    const std::array<int, 2> rb = { 0, 2 };
                    const scalar thick = 5.8889e-05 * 7.0;


                    for (const auto& col : rb) {
                        const scalar de_wp = 2.0 * M_PI / lambda[col] * thick * 1.55338;
                        const scalar do_ep = 2.0 * M_PI / lambda[col] * thick * 1.54425;

                        m(0, 0) = 0.5 * (exp(ii * de_wp) + exp(ii * do_ep));
                        m(0, 1) = 0.5 * (exp(ii * de_wp) - exp(ii * do_ep));
                        m(1, 0) = m(0, 1);
                        m(1, 1) = m(0, 0);

                        Eo[col] = m * Eo[col];
                    }

                }

                std::array<float, 4> color;

                for (int rgb = 0; rgb < 3; rgb++) {
                    scalar I = std::abs(Eo[rgb](0) * cos(th) + Eo[rgb](1) * sin(th));
                    I *= I;
                    color[rgb] = pow(intensity[rgb] * I, gamma);
                }

                color[3] = alpha;

                colorFunc(CData, color, cross2ind(i, j, (std::size_t)cumVoxProdSlice[0]));
            }
        }


}


}}}