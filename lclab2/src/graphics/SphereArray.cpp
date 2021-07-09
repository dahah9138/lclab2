#include "SphereArray.h"
#include <algorithm>

using namespace Magnum;

namespace LC {
    
void SphereArray::Init(void* positions, const Vector3& (*Access)(void* data, std::size_t i), std::size_t size, int subdivisions) {
    numObjects = size;
    sphereRadius = 0.5f / (Float)pow(numObjects, 1.0f/3.0f);
    spherePositions = Containers::Array<Vector3>{ NoInit, numObjects };
    sphereInstanceData = Containers::Array<SphereInstanceData>{ NoInit, numObjects };

    for (std::size_t i = 0; i < numObjects; ++i) {
        spherePositions[i] = Access(positions, i);

        /* Fill in the instance data. Most of this stays the same, except
           for the translation */
        sphereInstanceData[i].transformationMatrix =
            Matrix4::translation(spherePositions[i]) *
            Matrix4::scaling(Vector3{ sphereRadius });
        sphereInstanceData[i].normalMatrix =
            sphereInstanceData[i].transformationMatrix.normalMatrix();
        sphereInstanceData[i].color = Color3{1.0f, 1.0f, 1.0f};
    }

    sphereShader = Shaders::PhongGL{
                Shaders::PhongGL::Flag::VertexColor |
                Shaders::PhongGL::Flag::InstancedTransformation };
    sphereInstanceBuffer = GL::Buffer{};
    sphereMesh = MeshTools::compile(Primitives::icosphereSolid(subdivisions));
    sphereMesh.addVertexBufferInstanced(sphereInstanceBuffer, 1, 0,
        Shaders::PhongGL::TransformationMatrix{},
        Shaders::PhongGL::NormalMatrix{},
        Shaders::PhongGL::Color3{});
    sphereMesh.setInstanceCount(sphereInstanceData.size());
}

void SphereArray::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {
    using namespace Magnum;

    sphereInstanceBuffer.setData(sphereInstanceData, GL::BufferUsage::DynamicDraw);
    sphereShader
        .setProjectionMatrix(projection)
        .setTransformationMatrix(arcball->viewMatrix())
        .setNormalMatrix(arcball->viewMatrix().normalMatrix())
        .draw(sphereMesh);
}

}