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

    // Solver
    LC::Solver* _solver;
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


    _solver = new LC::FrankOseen::ElasticOnly::FOFDSolver;

    using Dataset = LC::FrankOseen::ElasticOnly::FOFDSolver::dataset;

    /* Setup data */
    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    data->voxels[0] = 10;
    data->voxels[1] = 10;
    data->voxels[2] = 10;

    data->cell_dims[0] = 1.0;
    data->cell_dims[2] = 1.0;
    data->cell_dims[3] = 1.0;

    data->k11 = LC::FrankOseen::ElasticConstants::_5CB("k11");
    data->k22 = LC::FrankOseen::ElasticConstants::_5CB("k22");
    data->k33 = LC::FrankOseen::ElasticConstants::_5CB("k33");

    _solver->Init();


    LC_INFO("Created sandbox!");
}

sandbox::~sandbox() {
	LC_INFO("Destroying sandbox.");

    delete _solver;
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


    // TODO: Add a button with imgui to start the relax
    bool pressedRelaxAndRelaxFinished = false;

    if (moving || pressedRelaxAndRelaxFinished) redraw();

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
