#ifndef SPHEREARRAY_H
#define SPHEREARRAY_H

#include "core.h"

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Trade/MeshData.h>
#include <Corrade/Containers/Optional.h>

#include "ArcBall.h"


namespace LC {

struct SphereArray {
	
    struct PolyInstanceData {
        Magnum::Matrix4 transformationMatrix;
        Magnum::Matrix3x3 normalMatrix;
        Magnum::Color3 color;
    };

    void Init(void* positions, std::function<Magnum::Vector3(void*, std::size_t)> Access, std::size_t size, int subdivisions = 2);
    void InitPositions(void* positions, const std::function<Magnum::Vector3(void*, std::size_t)> &Access, std::size_t size);
    void ZProfileColor(const float &hue1 = 240.f, const float &hue2 = 60.f);
    void Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4 &projection);
    void Draw(const Magnum::Matrix4& viewMatrix, const Magnum::Matrix4& projection);

    Magnum::UnsignedInt numObjects;

    Magnum::Containers::Array<Magnum::Vector3> polyPositions;
    Magnum::Float polyRadius = 1.0f;
    Magnum::GL::Mesh polyMesh{ Magnum::NoCreate };
    Magnum::GL::Buffer polyInstanceBuffer{ Magnum::NoCreate };
    Magnum::Shaders::PhongGL polyShader{ Magnum::NoCreate };
    Magnum::Containers::Array<PolyInstanceData> polyInstanceData;
};
}

#endif