#include "Application.h"

namespace LC
{
	Application::Application(const Arguments& arguments): Platform::Application{arguments} {
		if (!pfd::settings::available()) {
			LC_CORE_WARN("Portable File Dialogs are not available on this platform");
		}
	}
	Application::Application(const Arguments& arguments, const Configuration& configuration) : Platform::Application{ arguments, configuration } {
		if (!pfd::settings::available()) {
			LC_CORE_WARN("Portable File Dialogs are not available on this platform");
		}
	}

	void Application::setupCamera(const Float& param, CameraType cameraType) {
		using namespace Math::Literals;

		_cameraType = cameraType;

		if (cameraType == CameraType::Group) {
			_cameraObject
				.setParent(&_scene)
				.translate(Vector3::zAxis(param));
			(*(_camera = new SceneGraph::Camera3D{ _cameraObject }))
				.setAspectRatioPolicy(SceneGraph::AspectRatioPolicy::Extend)
				.setProjectionMatrix(Matrix4::perspectiveProjection(35.0_degf, 1.0f, 0.01f, 1000.0f))
				.setViewport(GL::defaultFramebuffer.viewport().size());

			/* Base object, parent of all (for easy manipulation) */
			_manipulator.setParent(&_scene);
		}
		else if (cameraType == CameraType::ArcBall) {
			const Vector3 eye = Vector3::zAxis(-5.0f);
			const Vector3 viewCenter;
			const Vector3 up = Vector3::yAxis();
			const Deg fov = 45.0_degf;
			_arcballCamera.emplace(eye, viewCenter, up, fov, windowSize());
			_arcballCamera->setLagging(param);

			_projectionMatrix = Matrix4::perspectiveProjection(fov,
				Vector2{ framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
		}
	}

	void Application::mouseScrollEvent(MouseScrollEvent& event) {
		_imgui.handleMouseScrollEvent(event);

		if (!_io->WantCaptureMouse) {


			if (_cameraType == CameraType::ArcBall) {
				const Float delta = event.offset().y();
				if (Math::abs(delta) >= 1.0e-2f) _arcballCamera->zoom(delta);
			}
			else if (_cameraType == CameraType::Group) {

				if (!event.offset().y()) return;

				/* Distance to origin */
				const Float distance = _cameraObject.transformation().translation().z();

				/* Move 15% of the distance back or forward */
				_cameraObject.translate(Vector3::zAxis(
					distance * (1.0f - (event.offset().y() > 0 ? 1 / 0.85f : 0.85f))));
			}
		}
		else _ioUpdate = true;

		event.setAccepted();
	}

	void Application::setupGUI() {



		/* Setup imgui (Calls ImGui::CreateContext() */
		_imgui = ImGuiIntegration::Context(Vector2{ windowSize() } / dpiScaling(),
			windowSize(), framebufferSize());

		// Create ImPlot Context() in pair with ImGui
		//if ((_options & App::OptionFlag::ImPlot) != App::OptionFlag::None)
		//	ImPlot::CreateContext();

		//LC::ImPlotIntegration::CreateContext();
		ImPlot::CreateContext();

		/* Set up proper blending to be used by ImGui. There's a great chance
		   you'll need this exact behavior for the rest of your scene. If not, set
		   this only for the drawFrame() call. */
		GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
			GL::Renderer::BlendEquation::Add);
		GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
			GL::Renderer::BlendFunction::OneMinusSourceAlpha);

		_io = &ImGui::GetIO();
	}

	void Application::viewportEvent(ViewportEvent& event) {
		GL::defaultFramebuffer.setViewport({ {}, event.framebufferSize() });

		if (_cameraType == CameraType::ArcBall) {
			_arcballCamera->reshape(event.windowSize());
			_projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
			Vector2{ event.framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
		}
		else if (_cameraType == CameraType::Group)
			_camera->setViewport(event.windowSize());

		_imgui.relayout(Vector2{ event.windowSize() } / event.dpiScaling(),
			event.windowSize(), event.framebufferSize());
	}

	void Application::mouseMoveEvent(MouseMoveEvent& event) {

		_imgui.handleMouseMoveEvent(event);

		if (!_io->WantCaptureMouse) {

			if (!event.buttons()) return;

			if (_cameraType == CameraType::Group) {

				if (event.buttons() & MouseMoveEvent::Button::Left) {
					const Vector3 currentPosition = positionOnSphere(event.position());
					const Vector3 axis = Math::cross(_previousPosition, currentPosition);

					if (_previousPosition.length() < 0.001f || axis.length() < 0.001f) return;

					_manipulator.rotate(Math::angle(_previousPosition, currentPosition), axis.normalized());
					_previousPosition = currentPosition;
				}
			}
			else if (_cameraType == CameraType::ArcBall) {
				if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
					_arcballCamera->translate(event.position());
				else _arcballCamera->rotate(event.position());
			}

		}
	}

	void Application::mousePressEvent(MouseEvent& event) {

		_imgui.handleMousePressEvent(event);

		if (!_io->WantCaptureMouse) {
			if (_cameraType == CameraType::ArcBall)
				_arcballCamera->initTransformation(event.position());
		}
		else _ioUpdate = true;
	}

	void Application::enableDepthTest() {
		GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
	}

	void Application::disableDepthTest() {
		GL::Renderer::disable(GL::Renderer::Feature::DepthTest);
	}

	void Application::enableFaceCulling() {
		GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
	}

	void Application::disableFaceCulling() {
		GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
	}

	void Application::guiRenderer() {
		GL::Renderer::enable(GL::Renderer::Feature::Blending);
		GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
		GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
		GL::Renderer::disable(GL::Renderer::Feature::DepthTest);
	}

	void Application::polyRenderer() {
		// Blending is needed for transparency
		//GL::Renderer::disable(GL::Renderer::Feature::Blending);
		//GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
		GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
		GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
	}

	Vector3 Application::positionOnSphere(const Vector2i& position) const {
		const Vector2 positionNormalized = Vector2{ position } / Vector2{ _camera->viewport() } -Vector2{ 0.5f };
		const Float length = positionNormalized.length();
		const Vector3 result(length > 1.0f ? Vector3(positionNormalized, 0.0f) : Vector3(positionNormalized, 1.0f - length));
		return (result * Vector3::yScale(-1.0f)).normalized();
	}

	Application::~Application() {

		// Good enough to just destroy ImPlot context here
		
		// Need to check that ImGui and ImPlot is being used
		bool condImGui = (_options & App::OptionFlag::ImGui) != App::OptionFlag::None;
		bool condImPlot = (_options & App::OptionFlag::ImPlot) != App::OptionFlag::None;

		if (condImGui && condImPlot) ImPlot::DestroyContext();

		Debug{} << "Terminating application.\n";
	}
}
