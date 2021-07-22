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
	
    struct SphereInstanceData {
        Magnum::Matrix4 transformationMatrix;
        Magnum::Matrix3x3 normalMatrix;
        Magnum::Color3 color;
    };

    void Init(void* positions, const Magnum::Vector3& (*Access)(void* data, std::size_t i), std::size_t size, int subdivisions = 2);
    void Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4 &projection);

    Magnum::UnsignedInt numObjects;

    Magnum::Containers::Array<Magnum::Vector3> spherePositions;
    Magnum::Float sphereRadius;
    Magnum::GL::Mesh sphereMesh{ Magnum::NoCreate };
    Magnum::GL::Buffer sphereInstanceBuffer{ Magnum::NoCreate };
    Magnum::Shaders::PhongGL sphereShader{ Magnum::NoCreate };
    Magnum::Containers::Array<SphereInstanceData> sphereInstanceData;
};
}

#endif