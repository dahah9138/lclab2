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
#include <Magnum/SceneGraph/Camera.h>

// ImGui
//#include <Magnum/ImGuiIntegration/Context.hpp>

// ImPlot addon for plotting
#include "implementation/ImContext.h"
//#include <implot.h>
#include "core.h"

#include "logger.h"
#include "graphics/ArcBall.h"
#include "solver/Solver.h"
#include "graphics/TransparentDrawable.h"

namespace LC
{
	using namespace Magnum;
	using std::chrono::high_resolution_clock;

	namespace App {
		enum class OptionFlag { None = 0, ImGui = BIT(0), ImPlot = BIT(1) };
		inline constexpr OptionFlag operator | (OptionFlag f1, OptionFlag f2) {
			return static_cast<OptionFlag>(static_cast<unsigned>(f1) | static_cast<unsigned>(f2));
		}
		inline constexpr OptionFlag operator & (OptionFlag f1, OptionFlag f2) {
			return static_cast<OptionFlag>(static_cast<unsigned>(f1) & static_cast<unsigned>(f2));
		}

	}

	class LC_API Application : public Platform::Application {
    public:

		// Not a flag
		enum class CameraType { ArcBall, Group };
		
        explicit Application(const Arguments& arguments);
		explicit Application(const Arguments& arguments, const Configuration& configuration);

		virtual void setupCamera(const Float& param, CameraType cameraType);
		virtual void mouseScrollEvent(MouseScrollEvent& event) override;
		virtual void viewportEvent(ViewportEvent& event) override;
		virtual void mouseMoveEvent(MouseMoveEvent& event) override;
		virtual void mousePressEvent(MouseEvent& event) override;

		void enableDepthTest();
		void disableDepthTest();

		void enableFaceCulling();
		void disableFaceCulling();

		Vector3 positionOnSphere(const Vector2i& position) const;

		void guiRenderer();
		void polyRenderer();

		void setupGUI();

		virtual ~Application();

		Containers::Optional<ArcBall> _arcballCamera;
		Matrix4 _projectionMatrix;

		// Make a part of Application
		Drawable::Scene3D _scene;
		Drawable::Object3D _manipulator, _cameraObject;
		SceneGraph::Camera3D* _camera;
		// Needed for camera
		Vector3 _previousPosition;

		CameraType _cameraType;
		App::OptionFlag _options = App::OptionFlag::ImGui;

		ImGuiIntegration::Context _imgui{ NoCreate };
		// Pointer to io context

		// Solver
    	std::unique_ptr<Solver> _solver;

		ImGuiIO* _io;
		bool _ioUpdate = true;

	};
	
}


#endif
