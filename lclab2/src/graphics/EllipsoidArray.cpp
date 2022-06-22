#include "EllipsoidArray.h"
#include <algorithm>

using namespace Magnum;

namespace LC {
    
void EllipsoidArray::Init(void* positions, std::function<Magnum::Vector3(void*, std::size_t)> Access, std::size_t size, int subdivisions) {
    numObjects = size;
    Float polyRadius = scale / (Float)pow(numObjects, 1.0f/3.0f);
    polyPositions = Containers::Array<Vector3>{ NoInit, numObjects };
    polyInstanceData = Containers::Array<PolyInstanceData>{ NoInit, numObjects };

    for (std::size_t i = 0; i < numObjects; ++i) {
        polyPositions[i] = Access(positions, i);

        /* Fill in the instance data. Most of this stays the same, except
           for the translation */
        polyInstanceData[i].transformationMatrix =
            Matrix4::translation(polyPositions[i]) * Matrix4::scaling(Vector3{ polyRadius });
        polyInstanceData[i].normalMatrix =
            polyInstanceData[i].transformationMatrix.normalMatrix();
        polyInstanceData[i].color = Color3{1.0f, 1.0f, 1.0f};
    }

    polyShader = Shaders::PhongGL{
                Shaders::PhongGL::Flag::VertexColor |
                Shaders::PhongGL::Flag::InstancedTransformation, 6 };
    polyInstanceBuffer = GL::Buffer{};

    polyMesh = MeshTools::compile(Primitives::capsule3DSolid(4, 4, 12, hLength, Magnum::Primitives::CapsuleTextureCoords()));
    polyMesh.addVertexBufferInstanced(polyInstanceBuffer, 1, 0,
        Shaders::PhongGL::TransformationMatrix{},
        Shaders::PhongGL::NormalMatrix{},
        Shaders::PhongGL::Color3{});
    polyMesh.setInstanceCount(polyInstanceData.size());
}

void EllipsoidArray::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {
    using namespace Magnum;

    polyInstanceBuffer.setData(polyInstanceData, GL::BufferUsage::DynamicDraw);
    polyShader
        .setSpecularColor({0.f, 0.f, 0.f, 0.f})
        .setProjectionMatrix(projection)
        .setTransformationMatrix(arcball->viewMatrix())
        .setNormalMatrix(arcball->viewMatrix().normalMatrix())
        .draw(polyMesh);
}

void EllipsoidArray::Draw(const Magnum::Matrix4& viewMatrix, const Magnum::Matrix4& projection) {
    using namespace Magnum;

    polyInstanceBuffer.setData(polyInstanceData, GL::BufferUsage::DynamicDraw);
    polyShader
        .setProjectionMatrix(projection)
        .setSpecularColor({ 0.f, 0.f, 0.f, 0.f })
        .setTransformationMatrix(viewMatrix)
        .setNormalMatrix(viewMatrix.normalMatrix())
        .draw(polyMesh);
}

}