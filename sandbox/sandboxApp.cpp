#include <lclab2.h>

using namespace Magnum;

/*
    Next: Maybe consider making a POM imaging lib for lclab2 (make with CPU first)
    - Needs Jones Matrix (I have matrices through Magnum)
    - Needs Lamp specifications
    - LC Optical properties
*/

class sandbox : public LC::Application
{
public:

	explicit sandbox(const Arguments& arguments);

	~sandbox();

private:
	virtual void drawEvent() override;

    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;

    // Tested geometries
    LC::SphereArray _grid;
    LC::Torus _sheet;
    LC::NormalTorus _sheetNormal;
};

sandbox::sandbox(const Arguments& arguments) : LC::Application{ arguments, Configuration{}.setTitle("Sandbox Application") } {

    /* Setup window and parameters */
    enableDepthTest();
    enableFaceCulling();

    /* Loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
    

    /* Setup camera */
    setupCamera(0.9f);

    /* Setup spheres */
    //_grid.Init();
    //_sheet.Init();
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

    //for (std::size_t i = 0; i < _grid.spherePositions.size(); ++i) {
    //    const Vector3 tmpCol = Vector3(std::rand(), std::rand(), std::rand()) /
    //        Float(RAND_MAX);
    //    _grid.sphereInstanceData[i].color = tmpCol;
    //}

    _arcballCamera->initTransformation(event.position());
    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void sandbox::mouseReleaseEvent(MouseEvent& event) {

    SDL_CaptureMouse(SDL_FALSE);
}


LC::Application* LC::createApplication(int argc, char **argv) {

	return new sandbox{ Platform::Application::Arguments{argc, argv} };
}
