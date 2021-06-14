#include "Application.h"

namespace LC
{
	Application::Application(const Arguments& arguments): Platform::Application{arguments} {
	

	}
	Application::Application(const Arguments& arguments, const Configuration& configuration) : Platform::Application{ arguments, configuration } {

	}

	void Application::setupCamera(const Float &lag) {

		using namespace Math::Literals;

		/* Setup camera */
		{
			const Vector3 eye = Vector3::zAxis(-5.0f);
			const Vector3 viewCenter;
			const Vector3 up = Vector3::yAxis();
			const Deg fov = 45.0_degf;
			_arcballCamera.emplace(eye, viewCenter, up, fov, windowSize());
			_arcballCamera->setLagging(lag);

			_projectionMatrix = Matrix4::perspectiveProjection(fov,
				Vector2{ framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
		}
	}

	void Application::mouseScrollEvent(MouseScrollEvent& event) {

		_imgui.handleMouseScrollEvent(event);

		if (!_io->WantCaptureMouse) {

			const Float delta = event.offset().y();

			if (Math::abs(delta) >= 1.0e-2f)
				_arcballCamera->zoom(delta);
		}
		else _ioUpdate = true;

		event.setAccepted();
	}

	void Application::viewportEvent(ViewportEvent& event) {
		GL::defaultFramebuffer.setViewport({ {}, event.framebufferSize() });
		_arcballCamera->reshape(event.windowSize());

		_projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
			Vector2{ event.framebufferSize() }.aspectRatio(), 0.01f, 100.0f);

		_imgui.relayout(Vector2{ event.windowSize() } / event.dpiScaling(),
			event.windowSize(), event.framebufferSize());
	}

	void Application::mouseMoveEvent(MouseMoveEvent& event) {

		if (!event.buttons()) return;

		_imgui.handleMouseMoveEvent(event);

		if (!_io->WantCaptureMouse) {
			if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
				_arcballCamera->translate(event.position());
			else _arcballCamera->rotate(event.position());
		}
		else _ioUpdate = true;

	}

	void Application::mousePressEvent(MouseEvent& event) {

		_imgui.handleMousePressEvent(event);

		if (!_io->WantCaptureMouse) {
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
		GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
		GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
		GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
		GL::Renderer::disable(GL::Renderer::Feature::Blending);
	}

	Application::~Application() {
		Debug{} << "Terminating application.\n";
	}
}
