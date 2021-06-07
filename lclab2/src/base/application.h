#ifndef LCLAB_H
#define LCLAB_H


#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/icoSphere.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Trade/MeshData.h>



#include "core.h"
#include "logger.h"

namespace LC
{
	using namespace Magnum;
	
	class LC_API application: public Platform::Application {
    public:
        explicit application(const Arguments& arguments);
		explicit application(const Arguments& arguments, const Configuration& configuration);

		virtual ~application();
	};
	
}


#endif