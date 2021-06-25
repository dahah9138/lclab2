#ifndef DYNAMIC_COLOR_SHEET_H
#define DYNAMIC_COLOR_SHEET_H

#include "core.h"
#include "logger.h"

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Corrade/Containers/Optional.h>

#include "ArcBall.h"

/*
	Polymorphic sheet class that will be inherited by any objects that can be made from a single sheet.
*/

namespace LC {

	struct LC_API DynamicColorSheet {

		typedef Magnum::Vector3(*PositionFunction)(Magnum::Float, Magnum::Float, Magnum::Float);

		struct Vertex {
			Magnum::Vector3 position;
			Magnum::Color4 color;
		};

		virtual void Init(PositionFunction pos = 0);
		virtual void Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection);

		// numpoints in x and y for sheet
		Magnum::UnsignedInt NX;
		Magnum::UnsignedInt NY;
		Magnum::Float CX;
		Magnum::Float CY;

		Magnum::GL::Mesh sheetMesh{ Magnum::NoCreate };
		Magnum::GL::Buffer sheetBuffer{ Magnum::NoCreate };
		Magnum::GL::Buffer sheetIndexBuffer{ Magnum::NoCreate };

		Magnum::Shaders::VertexColorGL3D sheetShader{ Magnum::NoCreate };
		Magnum::Containers::Array<Vertex> data;
		Magnum::Containers::Array<char> ind_data;

	};

}

#endif