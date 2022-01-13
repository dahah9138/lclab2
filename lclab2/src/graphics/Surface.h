#ifndef GRAPHICS_SURFACE_H
#define GRAPHICS_SURFACE_H

#include "core.h"
#include "logger.h"

#include <Magnum/ImageView.h>
#include <Magnum/Mesh.h>

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>

#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Compile.h>

#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>

#include <Magnum/Trade/MeshData.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Corrade/Containers/Optional.h>

#include "ArcBall.h"

/*
	Polymorphic sheet class that will be inherited by any objects that can be made from a single sheet.
*/

namespace LC {

    using namespace Magnum;

    struct Surface {

        struct Vertex {
            void operator = (const Vertex& vert) {
                position = vert.position;
                normal = vert.normal;
                color = vert.color;
            }
            Magnum::Vector3 position;
            Magnum::Vector3 normal;
            Magnum::Color4 color;
        };

        virtual void Init(Vertex* verts, unsigned int nVerts, unsigned int* inds, unsigned int nIndices, const Magnum::Vector3& translate = Magnum::Vector3{0.0f,0.0f,0.0f});

        Trade::MeshData Data();
        GL::Mesh Mesh();

        // This is the data to hold on to!
        Containers::Array<Vertex> vertices;
        Containers::Array<UnsignedInt> indices;
        Containers::Array<Trade::MeshAttributeData> attributes;

        // Hold onto vertex buffer
        GL::Buffer vertexBuffer;
    };

}

#endif