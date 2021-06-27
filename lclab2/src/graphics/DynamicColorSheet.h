#ifndef DYNAMIC_COLOR_SHEET_H
#define DYNAMIC_COLOR_SHEET_H

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
#include <Corrade/Containers/Optional.h>

#include "ArcBall.h"

/*
	Polymorphic sheet class that will be inherited by any objects that can be made from a single sheet.
*/

namespace LC {

    using namespace Magnum;

    struct LC_API DynamicColorSheet {

        typedef Magnum::Vector3(*PositionFunction)(Magnum::Float, Magnum::Float, Magnum::Float);

        struct Vertex {
            Magnum::Vector3 position;
            Magnum::Color4 color;
        };

        virtual void Init(PositionFunction pos = 0, Float offset = 0.0f);
        Trade::MeshData Data();
        GL::Mesh Mesh();

        // numpoints in x and y for sheet
        UnsignedInt NX = 32;
        UnsignedInt NY = 32;
        Float CX = 1.0f;
        Float CY = 1.0f;

        // This is the data to hold on to!
        Containers::Array<Vertex> vertices;
        Containers::Array<UnsignedInt> indices;
        Containers::Array<Trade::MeshAttributeData> attributes;


        // Hold onto vertex buffer
        GL::Buffer vertexBuffer;
    };

}

#endif