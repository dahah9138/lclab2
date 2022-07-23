#ifndef ZPROFILE_H
#define ZPROFILE_H

#include <lclab2.h>

using namespace Magnum;
using namespace Math::Literals;

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


#endif