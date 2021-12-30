#include "Application.h"

namespace LC
{
	Application::Application(const Arguments& arguments): Platform::Application{arguments} {

		_relaxFuture.second = false;

		if (!pfd::settings::available()) {
			LC_CORE_WARN("Portable File Dialogs are not available on this platform");
		}
	}
	Application::Application(const Arguments& arguments, const Configuration& configuration) : Platform::Application{ arguments, configuration } {
		
		_relaxFuture.second = false;
		
		if (!pfd::settings::available()) {
			LC_CORE_WARN("Portable File Dialogs are not available on this platform");
		}
	}

	void Application::setupCamera(const Float& param, CameraType cameraType) {
		using namespace Magnum::Math::Literals;

		_cameraType = cameraType;

		if (static_cast<int>(cameraType) & static_cast<int>(CameraType::Group)) {
			_cameraObject
				.setParent(&_scene)
				.translate(Vector3::zAxis(param));
			(*(_camera = new SceneGraph::Camera3D{ _cameraObject }))
				.setAspectRatioPolicy(SceneGraph::AspectRatioPolicy::Extend)
				.setProjectionMatrix(Matrix4::perspectiveProjection(35.0_degf, 1.0f, 0.01f, 1000.0f))
				.setViewport(GL::defaultFramebuffer.viewport().size());

			_manipulator = std::make_unique<Drawable::Object3D>();
			/* Base object, parent of all (for easy manipulation) */
			_manipulator->setParent(&_scene);
		}
		if (static_cast<int>(cameraType) & static_cast<int>(CameraType::ArcBall)) {
			const Vector3 eye = Vector3::zAxis(5.0f);
			const Vector3 viewCenter;
			const Vector3 up = Vector3::yAxis();
			const Deg fov = 45.0_degf;
			_arcballCamera.emplace(eye, viewCenter, up, fov, windowSize());
			_arcballCamera->setLagging(0.0);

			_projectionMatrix = Matrix4::perspectiveProjection(fov,
				Vector2{ framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
		}
	}

	void Application::mouseScrollEvent(MouseScrollEvent& event) {
		_imgui.handleMouseScrollEvent(event);

		if (!_io->WantCaptureMouse) {


			if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::ArcBall)) {
				const Float delta = event.offset().y();
				if (Magnum::Math::abs(delta) >= 1.0e-2f) _arcballCamera->zoom(delta);
			}
			if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::Group)) {

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
		if ((_options & App::OptionFlag::ImPlot) != App::OptionFlag::None)
			ImPlot::CreateContext();

		GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
			GL::Renderer::BlendEquation::Add);
		GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
			GL::Renderer::BlendFunction::OneMinusSourceAlpha);

		_io = &ImGui::GetIO();
	}

	bool Application::checkRelax(const std::size_t& seconds) {

		// Check if future ready
		if (_relaxFuture.second) {
			auto status = _relaxFuture.first.wait_for(std::chrono::seconds(seconds));

			if (status == std::future_status::ready) {
				_relaxFuture.first.get();
				_relaxFuture.second = false;
			}
		}
		return !_relaxFuture.second;
	}

	std::string Application::saveDialog() {
		auto sf = pfd::save_file("Select save file name", LCLAB2_ROOT_PATH,
			{ "LMT Files (.lmt .lmat)", "*.lmt *.lmat",
			  "All Files", "*" },
			pfd::opt::none);
		return sf.result();
	}

	std::vector<std::string> Application::openDialog() {
		auto of = pfd::open_file("Select save file name", LCLAB2_ROOT_PATH,
			{ "LMT Files (.lmt .lmat)", "*.lmt *.lmat",
					"All Files", "*" },
			pfd::opt::none);
		return of.result();
	}

	void Application::viewportEvent(ViewportEvent& event) {
		GL::defaultFramebuffer.setViewport({ {}, event.framebufferSize() });

		if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::ArcBall)) {
			_arcballCamera->reshape(event.windowSize());
			_projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
			Vector2{ event.framebufferSize() }.aspectRatio(), 0.01f, 100.0f);
		}
		if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::Group))
			_camera->setViewport(event.windowSize());

		_imgui.relayout(Vector2{ event.windowSize() } / event.dpiScaling(),
			event.windowSize(), event.framebufferSize());
	}

	void Application::mouseMoveEvent(MouseMoveEvent& event) {

		_imgui.handleMouseMoveEvent(event);

		if (!_io->WantCaptureMouse) {

			if (!event.buttons()) return;

			if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::Group)) {

				if (event.buttons() & MouseMoveEvent::Button::Left) {
					const Vector3 currentPosition = positionOnSphere(event.position());
					const Vector3 axis = Magnum::Math::cross(_previousPosition, currentPosition);

					if (_previousPosition.length() < 0.001f || axis.length() < 0.001f) return;

					_manipulator->rotate(Magnum::Math::angle(_previousPosition, currentPosition), axis.normalized());
					_previousPosition = currentPosition;
				}
			}
			if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::ArcBall)) {
				if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
					_arcballCamera->translate(event.position());
				else _arcballCamera->rotate(event.position());
			}

		}
	}

	void Application::mousePressEvent(MouseEvent& event) {

		if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::Group))
			if (event.button() == MouseEvent::Button::Left)
				_previousPosition = positionOnSphere(event.position());

		_imgui.handleMousePressEvent(event);

		if (static_cast<int>(_cameraType) & static_cast<int>(CameraType::ArcBall))
			if (!_io->WantCaptureMouse)
				_arcballCamera->initTransformation(event.position());
	}

	void Application::mouseReleaseEvent(MouseEvent& event) {
		if (event.button() == MouseEvent::Button::Left)
			_previousPosition = Vector3();
		if (_imgui.handleMouseReleaseEvent(event)) _ioUpdate = true;
	}

	void Application::keyPressEvent(KeyEvent& event) {
		// Check if Ctrl + S or Ctrl + O is pressed
		if ((event.key() == KeyEvent::Key::S) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { save(); }
		else if ((event.key() == KeyEvent::Key::O) && (event.modifiers() & KeyEvent::Modifier::Ctrl)) { open(); }
		if (_imgui.handleKeyPressEvent(event)) _ioUpdate = true;
	}

	void Application::keyReleaseEvent(KeyEvent& event) {
		if ((event.key() == KeyEvent::Key::S)) {}
		else if ((event.key() == KeyEvent::Key::O)) {}
		if (_imgui.handleKeyReleaseEvent(event)) _ioUpdate = true;
	}

	void Application::textInputEvent(TextInputEvent& event) {
		_imgui.handleTextInputEvent(event);
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

	void Application::sortObjects(Magnum::SceneGraph::DrawableGroup3D& drawables) {
		// Sort objects to draw in correct order
		std::vector<std::pair<std::reference_wrapper<SceneGraph::Drawable3D>, Matrix4>>
			drawableTransformations = _camera->drawableTransformations(drawables);

		std::sort(drawableTransformations.begin(), drawableTransformations.end(),
			[](const std::pair<std::reference_wrapper<SceneGraph::Drawable3D>, Matrix4>& a,
				const std::pair<std::reference_wrapper<SceneGraph::Drawable3D>, Matrix4>& b) {
					return a.second.translation().z() < b.second.translation().z();
			});
	}

	void Application::saveMenu(bool& loaded, std::function<void()> loadAction) {
		if (ImGui::BeginMenuBar()) {
			if (ImGui::BeginMenu("File")) {

				//if (ImGui::MenuItem("New", "Ctrl+N")) {}
				if (ImGui::MenuItem("Open", "Ctrl+O")) {
					loaded = open();
					if (loaded && loadAction) loadAction();
					else if (loaded) LC_CORE_WARN("menu load action was null");
				}
				if (ImGui::MenuItem("Save", "Ctrl+S")) {
					save();
				}
				if (ImGui::MenuItem("Save As..")) {
					saveAs();
				}

				ImGui::EndMenu();
			}
			ImGui::EndMenuBar();
		}
	}

	void Application::save() {
		bool ready = checkRelax();
		if (!ready) return;

		// Check if file exists (was loaded)

		std::string res;

		bool loadedFromFile = !_header.readFile.empty() || !_header.writeFile.empty();

		if (!loadedFromFile) {
			res = saveDialog();
		}
		else {
			res = _header.readFile;
		}

		if (!res.empty()) {
			_header.writeFile = res;
			_solver->Export(_header);
			LC_CORE_INFO("Saved to file {0}", res.c_str());
			_header.readFile = res;
		}
	}
	void Application::saveAs() {
		bool ready = checkRelax();
		if (!ready) return;

		// Pass file to header to be written
		auto res = saveDialog();

		if (!res.empty()) {
			_header.writeFile = res;
			_solver->Export(_header);
			LC_INFO("Saved to file {0}", res.c_str());
			_header.readFile = res;
		}
	}
	bool Application::open() {
		bool ready = checkRelax();
		if (!ready) return false;
		bool success = false;

		auto res = openDialog();

		if (!res.empty()) {
			_header.readFile = res[0];
			_solver->Import(_header);
			success = true;
			LC_CORE_INFO("Loaded file {0}", res[0].c_str());

		}

		return success;
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
