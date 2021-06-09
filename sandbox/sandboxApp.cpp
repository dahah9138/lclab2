#include <lclab2.h>

using namespace Magnum;

/*
    Next: Maybe consider making a POM imaging lib for lclab2 (make with CPU first)
*/

class sandbox : public LC::application
{
public:

	explicit sandbox(const Arguments& arguments);

	~sandbox();

private:
	virtual void drawEvent() override;

    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;
    void mouseMoveEvent(MouseMoveEvent& event) override;
    void viewportEvent(ViewportEvent& event) override;
    void mouseScrollEvent(MouseScrollEvent& event) override;

    Containers::Optional<ArcBall> _arcballCamera;

    LC::SphereArray _grid;
    LC::Torus _sheet;
    LC::NormalTorus _sheetNormal;

    Matrix4 _projectionMatrix;

};

sandbox::sandbox(const Arguments& arguments) : LC::application{ arguments, Configuration{}.setTitle("Sandbox Application") } {
	
    using namespace Math::Literals;

    /* Setup window and parameters */
    {
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

        /* Loop at 60 Hz max */
        setSwapInterval(1);
        setMinimalLoopPeriod(16);
    }

    /* Setup camera */
    {
        const Vector3 eye = Vector3::zAxis(5.0f);
        const Vector3 viewCenter;
        const Vector3 up = Vector3::yAxis();
        const Deg fov = 45.0_degf;
        _arcballCamera.emplace(eye, viewCenter, up, fov, windowSize());
        _arcballCamera->setLagging(0.55f);

        _projectionMatrix = Matrix4::perspectiveProjection(fov,
            Vector2{ framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
    }

    /* Setup spheres */
    _grid.Init();
    _sheet.Init();
    _sheetNormal.Init();

    LC_INFO("Created sandbox!");

}

sandbox::~sandbox() {
	LC_INFO("Destroying sandbox.");
}

void sandbox::drawEvent()
{
    GL::defaultFramebuffer.clear(
        GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    /* Update camera before drawing instances */
    const bool moving = _arcballCamera->updateTransformation();

    //_grid.Draw(_arcballCamera, _projectionMatrix);
    //_sheet.Draw(_arcballCamera, _projectionMatrix);
    _sheetNormal.Draw(_arcballCamera, _projectionMatrix);

    swapBuffers();

    if (moving) redraw();
}

void sandbox::mousePressEvent(MouseEvent& event) {

    // Configure the mouse press event
    SDL_CaptureMouse(SDL_TRUE);

    for (std::size_t i = 0; i < _grid.spherePositions.size(); ++i) {
        const Vector3 tmpCol = Vector3(std::rand(), std::rand(), std::rand()) /
            Float(RAND_MAX);
        _grid.sphereInstanceData[i].color = tmpCol;
    }


    _arcballCamera->initTransformation(event.position());
    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void sandbox::mouseReleaseEvent(MouseEvent& event) {

    SDL_CaptureMouse(SDL_FALSE);
}

void sandbox::mouseMoveEvent(MouseMoveEvent& event) {

    if (!event.buttons()) return;

    if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
        _arcballCamera->translate(event.position());
    else _arcballCamera->rotate(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void sandbox::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({ {}, event.framebufferSize() });
    _arcballCamera->reshape(event.windowSize());

    _projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
        Vector2{ event.framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
}

void sandbox::mouseScrollEvent(MouseScrollEvent& event) {
    const Float delta = event.offset().y();
    if (Math::abs(delta) < 1.0e-2f) return;

    _arcballCamera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

LC::application* LC::createApplication(int argc, char **argv)
{

	return new sandbox{ Platform::Application::Arguments{argc, argv} };
}
