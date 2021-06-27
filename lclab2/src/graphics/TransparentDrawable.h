#ifndef TRANSPARENT_DRAWABLE_H
#define TRANSPARENT_DRAWABLE_H

#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>

#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>

#include <Magnum/Shaders/VertexColor.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>

namespace LC {	namespace Drawable {

	using namespace Magnum;

	typedef SceneGraph::Object<SceneGraph::MatrixTransformation3D> Object3D;
	typedef SceneGraph::Scene<SceneGraph::MatrixTransformation3D> Scene3D;

    class TransparentFlatDrawable : public SceneGraph::Drawable3D {
    public:
        explicit TransparentFlatDrawable(Object3D& object, Shaders::VertexColor3D& shader, GL::Mesh& mesh, const Vector3& position, SceneGraph::DrawableGroup3D& group) : SceneGraph::Drawable3D{ object, &group }, _shader(shader), _mesh(mesh), _position{ position } {}

    private:
        void draw(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera) override;

        Shaders::VertexColorGL3D& _shader;
        GL::Mesh& _mesh;
        Vector3 _position;
    };

    // Really a draw call for FlatTransparentDrawable
    void TransparentFlatDrawable::draw(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera) {
        GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::enable(GL::Renderer::Feature::Blending);
        _shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix * Matrix4::translation(_position))
            .draw(_mesh);
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::disable(GL::Renderer::Feature::Blending);
    }

}}


#endif