#ifndef ZPROFILE_H
#define ZPROFILE_H

#include <lclab2.h>

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

    void GenerateProfile(const LC::scalar *rz, std::array<int,3> vox, std::array<LC::scalar,3> cell) {
                
        graph = Graph{};

        graph->vox_x = vox[0];
        graph->vox_y = vox[1];

        unsigned int plane = vox[0] * vox[1];

        graph->data = std::unique_ptr<Magnum::Vector3[]>(new Magnum::Vector3[plane]);


        float dx = cell[0] / (vox[0] - 1);
        float dy = cell[1] / (vox[1] - 1);

        // Determine how large the spheres should be:
        graph->grid.polyRadius = 1.f/3.f * sqrt(dx * dx + dy * dy);

        // Find zmin and set as zero point
        LC::scalar zmin = rz[0];
        for (int i = 1; i < plane; i++) {
            if (zmin > rz[i])
                zmin = rz[i];
        }

        // Evaulate each point in the xy plane
        for (int x = 0; x < vox[0]; x++) {
            for (int y = 0; y < vox[1]; y++) {
                unsigned int idx = x + y * vox[0];
                // Initialize data
                graph->data[idx] = Magnum::Vector3( -cell[0]/2. + x * dx, -cell[1] / 2. + y * dy, rz[idx] - zmin);
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
    void Draw(const Containers::Optional<LC::ArcBall>& arcball, const Magnum::Matrix4& projectionMatrix) {
        graph->grid.Draw(arcball, projectionMatrix);
    }
    Containers::Optional<Graph> graph;
    Box box;
};

#endif