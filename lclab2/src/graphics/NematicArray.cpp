#include "NematicArray.h"
#include <Magnum/GL/Renderer.h>
#include <algorithm>

using namespace Magnum;

namespace LC {
    
void NematicArray::Init(void* positions, std::function<Magnum::Vector3(void*, std::size_t)> Access, std::size_t size, int subdivisions) {
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
        // Compute the normal matrix from the transformation matrix in advance
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

    // Initialize meshes
    polyMeshes.resize(map.size());

    polyMeshes[DrawType::Ellipsoid] = MeshTools::compile(Primitives::capsule3DSolid(4, 4, 12, hLength, Magnum::Primitives::CapsuleTextureCoords()));
    polyMeshes[DrawType::Cone] = MeshTools::compile(Primitives::coneSolid(3, 12, hLength, Primitives::ConeFlag::CapEnd));
    polyMeshes[DrawType::Cylinder] = MeshTools::compile(Primitives::cylinderSolid(3, 12, hLength, Primitives::CylinderFlag::CapEnds));

    // Instantiate meshes
    for (int i = 0; i < map.size(); i++) {
        polyMeshes[i].addVertexBufferInstanced(polyInstanceBuffer, 1, 0,
            Shaders::PhongGL::TransformationMatrix{},
            Shaders::PhongGL::NormalMatrix{},
            Shaders::PhongGL::Color3{});
        polyMeshes[i].setInstanceCount(polyInstanceData.size());
    }
}

void NematicArray::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {
    using namespace Magnum;
    GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
    polyInstanceBuffer.setData(polyInstanceData, GL::BufferUsage::DynamicDraw);
    polyShader
        .setSpecularColor({0.f, 0.f, 0.f, 0.f})
        .setProjectionMatrix(projection)
        .setTransformationMatrix(arcball->viewMatrix())
        .setNormalMatrix(arcball->viewMatrix().normalMatrix())
        .draw(polyMeshes[selected_drawType]);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
}

void NematicArray::Draw(const Magnum::Matrix4& viewMatrix, const Magnum::Matrix4& projection) {
    using namespace Magnum;
    GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
    polyInstanceBuffer.setData(polyInstanceData, GL::BufferUsage::DynamicDraw);
    polyShader
        .setProjectionMatrix(projection)
        .setSpecularColor({ specular, specular, specular, specular })
        .setDiffuseColor({ diffuse, diffuse, diffuse, diffuse })
        .setAmbientColor({ ambient, ambient, ambient, ambient })
        .setTransformationMatrix(viewMatrix)
        .setNormalMatrix(viewMatrix.normalMatrix())
        .draw(polyMeshes[selected_drawType]);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
}

}