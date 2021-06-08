#include <lclab2.h>

using namespace Magnum;

/*
    Currently just Magnum tutorial 2.
    TODO:
    1. Create a surface primitive
    2. Make a struct called Vertex with position and color
    3. Feed this container of vertices to the surface primitive
    4. Use the shader settings to color blend the plane
    5. Port plane primitive (with Vertex declared inside class) to lclab2 library
    6. Call from lclab2
    7. Repeat until all desired primitives have been created (e.g. cylinder, arrow, ...)
    
    Next: Maybe consider making a POM imaging lib for lclab2 (make with CPU first)
*/

struct SphereInstanceData {
    Matrix4 transformationMatrix;
    Matrix3x3 normalMatrix;
    Color3 color;
};

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

	GL::Mesh _mesh;
    Shaders::PhongGL _shader;
    Containers::Optional<ArcBall> _arcballCamera;

    Containers::Array<Vector3> _spherePositions;
    Float _sphereRadius;

    GL::Mesh _sphereMesh{ NoCreate };
    GL::Buffer _sphereInstanceBuffer{ NoCreate };
    Shaders::PhongGL _sphereShader{ NoCreate };
    Containers::Array<SphereInstanceData> _sphereInstanceData;


	Matrix4 _projectionMatrix;

};

sandbox::sandbox(const Arguments& arguments) : LC::application{ arguments, Configuration{}.setTitle("Cube test") } {
	
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
    {
        _sphereRadius = 0.03f;
        const UnsignedInt NX = 10;
        const UnsignedInt NY = 10;
        const UnsignedInt numSpheres = NX * NY;
        _spherePositions = Containers::Array<Vector3>{ NoInit, numSpheres };
        _sphereInstanceData = Containers::Array<SphereInstanceData>{ NoInit, numSpheres };

        for (std::size_t i = 0; i < numSpheres; ++i) {

            UnsignedInt ii = i / NY;
            UnsignedInt jj = i - ii * NY;
            Float x = (Float)ii / NX;
            Float y = (Float)jj / NY;

            _spherePositions[i] = Vector3(-0.5f + x, -0.5f + y, 0.0f);

            //_spherePositions[i] = tmpPos * 2.0f - Vector3{ 1.0f };
            //_spherePositions[i].y() *= 0.5f;

            /* Fill in the instance data. Most of this stays the same, except
               for the translation */
            _sphereInstanceData[i].transformationMatrix =
                Matrix4::translation(_spherePositions[i]) *
                Matrix4::scaling(Vector3{ _sphereRadius });
            _sphereInstanceData[i].normalMatrix =
                _sphereInstanceData[i].transformationMatrix.normalMatrix();
            _sphereInstanceData[i].color = Color3::green();
        }

        _sphereShader = Shaders::PhongGL{
                    Shaders::PhongGL::Flag::VertexColor |
                    Shaders::PhongGL::Flag::InstancedTransformation };
        _sphereInstanceBuffer = GL::Buffer{};
        _sphereMesh = MeshTools::compile(Primitives::icosphereSolid(2));
        _sphereMesh.addVertexBufferInstanced(_sphereInstanceBuffer, 1, 0,
            Shaders::PhongGL::TransformationMatrix{},
            Shaders::PhongGL::NormalMatrix{},
            Shaders::PhongGL::Color3{});
        _sphereMesh.setInstanceCount(_sphereInstanceData.size());
    }

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

    for (std::size_t i = 0; i != _spherePositions.size(); ++i)
        _sphereInstanceData[i].transformationMatrix.translation() =
        _spherePositions[i];

    _sphereInstanceBuffer.setData(_sphereInstanceData, GL::BufferUsage::DynamicDraw);
    _sphereShader
        .setProjectionMatrix(_projectionMatrix)
        .setTransformationMatrix(_arcballCamera->viewMatrix())
        .setNormalMatrix(_arcballCamera->viewMatrix().normalMatrix())
        .draw(_sphereMesh);

    swapBuffers();

    if (moving) redraw();
}

void sandbox::mousePressEvent(MouseEvent& event) {

    // Configure the mouse press event
    SDL_CaptureMouse(SDL_TRUE);

    for (std::size_t i = 0; i < _spherePositions.size(); ++i) {
        const Vector3 tmpCol = Vector3(std::rand(), std::rand(), std::rand()) /
            Float(RAND_MAX);
        _sphereInstanceData[i].color = tmpCol;
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
