#ifndef APPLICATION_H
#define APPLICATION_H


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
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Trade/MeshData.h>



#include "core.h"
#include "logger.h"
#include "../graphics/ArcBall.h"

namespace LC
{
	using namespace Magnum;
	
	class LC_API Application : public Platform::Application {
    public:
        explicit Application(const Arguments& arguments);
		explicit Application(const Arguments& arguments, const Configuration& configuration);

		virtual void setupCamera(const Float& lag);
		virtual void mouseScrollEvent(MouseScrollEvent& event) override;
		virtual void viewportEvent(ViewportEvent& event) override;
		virtual void mouseMoveEvent(MouseMoveEvent& event) override;
		virtual void mousePressEvent(MouseEvent& event) override;

		inline void enableDepthTest();
		inline void disableDepthTest();

		inline void enableFaceCulling();
		inline void disableFaceCulling();

		virtual ~Application();

		Containers::Optional<ArcBall> _arcballCamera;
		Matrix4 _projectionMatrix;
	};
	
}


#endif
