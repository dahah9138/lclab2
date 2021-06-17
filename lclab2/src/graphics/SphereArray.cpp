#include "SphereArray.h"
#include <algorithm>

namespace LC {
void SphereArray::Init() {
    using namespace Magnum;

    sphereRadius = 0.5f/(Float)((std::max)(NX, NY)-1);
    UnsignedInt numSpheres = NX * NY;
    spherePositions = Containers::Array<Vector3>{ NoInit, numSpheres };
    sphereInstanceData = Containers::Array<SphereInstanceData>{ NoInit, numSpheres };

    for (std::size_t i = 0; i < numSpheres; ++i) {

        // Use indexing matlab format
        UnsignedInt jj = i / NY;
        UnsignedInt ii = i - jj * NY;
        Float x = (Float)ii / (NX-1);
        Float y = (Float)jj / (NY-1);

        spherePositions[i] = Vector3(-0.5f + x, CY/CX*(-0.5f + y), 0.0f);

        /* Fill in the instance data. Most of this stays the same, except
           for the translation */
        sphereInstanceData[i].transformationMatrix =
            Matrix4::translation(spherePositions[i]) *
            Matrix4::scaling(Vector3{ sphereRadius });
        sphereInstanceData[i].normalMatrix =
            sphereInstanceData[i].transformationMatrix.normalMatrix();
        sphereInstanceData[i].color = Color3::green();
    }

    sphereShader = Shaders::PhongGL{
                Shaders::PhongGL::Flag::VertexColor |
                Shaders::PhongGL::Flag::InstancedTransformation };
    sphereInstanceBuffer = GL::Buffer{};
    sphereMesh = MeshTools::compile(Primitives::icosphereSolid(2));
    sphereMesh.addVertexBufferInstanced(sphereInstanceBuffer, 1, 0,
        Shaders::PhongGL::TransformationMatrix{},
        Shaders::PhongGL::NormalMatrix{},
        Shaders::PhongGL::Color3{});
    sphereMesh.setInstanceCount(sphereInstanceData.size());
}

void SphereArray::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {
    using namespace Magnum;

    for (std::size_t i = 0; i != spherePositions.size(); ++i)
        sphereInstanceData[i].transformationMatrix.translation() =
        spherePositions[i];

    sphereInstanceBuffer.setData(sphereInstanceData, GL::BufferUsage::DynamicDraw);
    sphereShader
        .setProjectionMatrix(projection)
        .setTransformationMatrix(arcball->viewMatrix())
        .setNormalMatrix(arcball->viewMatrix().normalMatrix())
        .draw(sphereMesh);
}

}