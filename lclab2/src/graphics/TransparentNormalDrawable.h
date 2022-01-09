#ifndef TRANSPARENT_NORMAL_DRAWABLE_H
#define TRANSPARENT_NORMAL_DRAWABLE_H

#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>

#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>

#include <Magnum/Shaders/PhongGL.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>

namespace LC {	namespace Drawable {

	using namespace Magnum;

	typedef SceneGraph::Object<SceneGraph::MatrixTransformation3D> Object3D;
	typedef SceneGraph::Scene<SceneGraph::MatrixTransformation3D> Scene3D;

    class TransparentNormalDrawable : public SceneGraph::Drawable3D {
    public:
        explicit TransparentNormalDrawable(Object3D& object, Shaders::PhongGL& shader, GL::Mesh& mesh, bool& draw, SceneGraph::DrawableGroup3D& group) : SceneGraph::Drawable3D{ object, &group }, _shader(shader), _mesh(mesh), _draw{draw} {}

    private:
        void draw(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera) override;

        Shaders::PhongGL& _shader;
        bool& _draw;
        GL::Mesh& _mesh;
    };

    // Really a draw call for FlatTransparentDrawable
    void TransparentNormalDrawable::draw(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera) {

        if (!_draw) return;

        //GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        //GL::Renderer::enable(GL::Renderer::Feature::Blending);
		_shader.setTransformationMatrix(transformationMatrix)
			   .setNormalMatrix(transformationMatrix.normalMatrix())
			   .setProjectionMatrix(camera.projectionMatrix())
			   .draw(_mesh);
        //GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        //GL::Renderer::disable(GL::Renderer::Feature::Blending);
    }

}}


#endif