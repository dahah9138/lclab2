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
        explicit TransparentNormalDrawable(Object3D& object, Shaders::PhongGL& shader, GL::Mesh& mesh, bool& draw, SceneGraph::DrawableGroup3D& group)
            : SceneGraph::Drawable3D{ object, &group }, _shader(shader), _mesh(mesh), _draw{draw} {}
        explicit TransparentNormalDrawable(Object3D& object, Shaders::PhongGL& shader, GL::Mesh& mesh, bool& draw, bool& noCull, SceneGraph::DrawableGroup3D& group)
            : SceneGraph::Drawable3D{ object, &group }, _shader(shader), _mesh(mesh), _draw{ draw }, _noFaceCull( &noCull ) {}
    private:
        void draw(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera) override;

        Shaders::PhongGL& _shader;
        bool* _noFaceCull = 0;
        bool& _draw;
        GL::Mesh& _mesh;
    };

    // Really a draw call for FlatTransparentDrawable
    void TransparentNormalDrawable::draw(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera) {

        if (!_draw) return;

        bool noCull = false;

        if (_noFaceCull) /* Check if exists */
            noCull = *_noFaceCull;

        GL::Renderer::enable(GL::Renderer::Feature::Multisampling);

        if (noCull)
            GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
        //GL::Renderer::enable(GL::Renderer::Feature::Blending);
        _shader.setTransformationMatrix(transformationMatrix)
               .setSpecularColor({0.f, 0.f, 0.f, 0.f})
			   .setNormalMatrix(transformationMatrix.normalMatrix())
			   .setProjectionMatrix(camera.projectionMatrix())
			   .draw(_mesh);
        if (noCull)
            GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

        GL::Renderer::disable(GL::Renderer::Feature::Multisampling);
        
        
    }

}}


#endif