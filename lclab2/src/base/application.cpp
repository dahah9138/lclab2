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
		const Float delta = event.offset().y();
		if (Math::abs(delta) < 1.0e-2f) return;

		_arcballCamera->zoom(delta);

		event.setAccepted();
		redraw(); /* camera has changed, redraw! */

		LC_INFO("Overrided once");
	}

	void Application::viewportEvent(ViewportEvent& event) {
		GL::defaultFramebuffer.setViewport({ {}, event.framebufferSize() });
		_arcballCamera->reshape(event.windowSize());

		_projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
			Vector2{ event.framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
	}

	Application::~Application() {
		Debug{} << "Terminating application.\n";
	}
}