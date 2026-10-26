#ifndef MESH_PROCESSOR_H
#define MESH_PROCESSOR_H

#include <vector>
#include <array>
#include <memory>

#include "Isosurface.h"

namespace LC { namespace Math {
    class MeshProcessor {
    public:
        struct ExtractedMesh {
            std::vector<IsoVertex> vertices;
            std::vector<unsigned int> indices;

            bool empty() const { return vertices.empty() || indices.empty(); }
            size_t numVertices() const { return vertices.size(); }
            size_t numTriangles() const { return indices.size() / 3; }
        };

        struct ExtractedLines {
            std::vector<float> vertices;
            std::vector<unsigned int> lineIndices;

            bool empty() const { return vertices.empty() || lineIndices.empty(); }
            size_t numVertices() const { return vertices.size() / 3; }
            size_t numSegments() const { return lineIndices.size() / 2; }
        };

        MeshProcessor();
        MeshProcessor(int Nx, int Ny, int Nz, float isoValue = 0.0f);
        ~MeshProcessor();

        MeshProcessor(const MeshProcessor&) = delete;
        MeshProcessor& operator=(const MeshProcessor&) = delete;
        MeshProcessor(MeshProcessor&&) noexcept;
        MeshProcessor& operator=(MeshProcessor&&) noexcept;

        void SetResolution(int Nx, int Ny, int Nz);
        void SetSpacing(float dx, float dy, float dz);
        void SetOrigin(float ox, float oy, float oz);
        void SetIsoValue(float isoValue);

        ExtractedMesh ComputeIsosurface(const float* voxelData);
        ExtractedMesh ComputeIsosurface(const float* voxelData, float isoValue);

        ExtractedLines ComputeBasisPreimageLines(
            const float* directorField,
            const std::array<float, 3>& n0,
            float alignmentThreshold = 0.0f
        );

        ExtractedMesh ComputeBasisPreimageTube(
            const float* directorField,
            const std::array<float, 3>& n0,
            float tubeRadius = 0.5f,
            int tubeSides = 8,
            float alignmentThreshold = 0.0f
        );

        ExtractedLines ComputeUmbilicLines(
            const float* directorField,
            float alignmentThreshold = 0.2f
        );

        ExtractedMesh ComputeUmbilicTube(
            const float* directorField,
            float tubeRadius,
            int tubeSides = 12,
            float alignmentThreshold = 0.2f
        );

    private:
        struct Pipeline;
        std::unique_ptr<Pipeline> vtk;

        int nx = 0;
        int ny = 0;
        int nz = 0;
        float iso = 0.0f;

        static void BuildOrthonormalBasis(
            const std::array<float, 3>& n0,
            std::array<float, 3>& b1,
            std::array<float, 3>& b2
        );

        float SampleAlignmentTrilinear(
            const float* directorField,
            const std::array<float, 3>& n0,
            float x, float y, float z
        ) const;

        ExtractedMesh ExtractPolyData(void* polyDataPtr, bool flipOrientation) const;
        ExtractedMesh BuildTubeFromExtractedLines(const ExtractedLines& lines, float tubeRadius, int tubeSides);
    };

}}

#endif