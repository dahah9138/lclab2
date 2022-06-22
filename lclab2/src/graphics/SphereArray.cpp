#include "SphereArray.h"
#include <algorithm>

using namespace Magnum;

namespace LC {
    
void SphereArray::InitPositions(void* positions, const std::function<Magnum::Vector3(void*, std::size_t)> &Access, std::size_t size) {
    polyPositions = Containers::Array<Vector3>{ NoInit, numObjects };
    for (std::size_t i = 0; i < numObjects; ++i) {
        polyPositions[i] = Access(positions, i);

        /* Fill in the instance data. Most of this stays the same, except
           for the translation */
        polyInstanceData[i].transformationMatrix =
            Matrix4::translation(polyPositions[i]) * Matrix4::scaling(Vector3{ polyRadius });
        polyInstanceData[i].normalMatrix =
            polyInstanceData[i].transformationMatrix.normalMatrix();
        polyInstanceData[i].color = Color3{ 1.0f, 1.0f, 1.0f };
    }
}

void SphereArray::ZProfileColor(const float& hue1, const float& hue2) {
    // Search for zmin and zmax
    float zmin(polyPositions[0].z()), zmax(polyPositions[0].z());
    for (std::size_t i = 1; i < numObjects; ++i) {
        float val = polyPositions[i].z();

        zmin = zmin > val ? val : zmin;
        zmax = zmax < val ? val : zmax;
    }


    // Linearly color according to scheme
    for (std::size_t i = 0; i < numObjects; ++i) {
        // shift by zmin and normalize z by zmax + zmin

        float zbar;
        if (zmax != zmin)
            zbar = (polyPositions[i].z() - zmin) / (zmax - zmin);
        else
            zbar = 1.0;

        float hue = (1.0f - zbar) * hue1 + zbar * hue2;

        polyInstanceData[i].color = Magnum::Color3::fromHsv(Deg{ hue }, 1.0f, 1.0f);
        // Normalize color
        for (int d = 0; d < 3; d++)
            polyInstanceData[i].color[d] = polyInstanceData[i].color[d] > 1. ? 1. : polyInstanceData[i].color[d];
    }
}

void SphereArray::Init(void* positions, std::function<Magnum::Vector3(void*, std::size_t)> Access, std::size_t size, int subdivisions) {
    numObjects = size;
    polyInstanceData = Containers::Array<PolyInstanceData>{ NoInit, numObjects };

    InitPositions(positions, Access, size);
    ZProfileColor();

    polyShader = Shaders::PhongGL{
                Shaders::PhongGL::Flag::VertexColor |
                Shaders::PhongGL::Flag::InstancedTransformation };
    polyInstanceBuffer = GL::Buffer{};
    polyMesh = MeshTools::compile(Primitives::icosphereSolid(subdivisions));
    polyMesh.addVertexBufferInstanced(polyInstanceBuffer, 1, 0,
        Shaders::PhongGL::TransformationMatrix{},
        Shaders::PhongGL::NormalMatrix{},
        Shaders::PhongGL::Color3{});
    polyMesh.setInstanceCount(polyInstanceData.size());
}

void SphereArray::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {
    using namespace Magnum;

    polyInstanceBuffer.setData(polyInstanceData, GL::BufferUsage::DynamicDraw);
    polyShader
        .setProjectionMatrix(projection)
        .setTransformationMatrix(arcball->viewMatrix())
        .setNormalMatrix(arcball->viewMatrix().normalMatrix())
        .draw(polyMesh);
}

void SphereArray::Draw(const Magnum::Matrix4 &viewMatrix, const Magnum::Matrix4& projection) {
    using namespace Magnum;

    polyInstanceBuffer.setData(polyInstanceData, GL::BufferUsage::DynamicDraw);
    polyShader
        .setProjectionMatrix(projection)
        .setTransformationMatrix(viewMatrix)
        .setNormalMatrix(viewMatrix.normalMatrix())
        .draw(polyMesh);
}

}