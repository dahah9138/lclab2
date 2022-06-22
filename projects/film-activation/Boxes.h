#ifndef BOXES_H
#define BOXES_H

#include <lclab2.h>

//Trade::MeshData primitiveVolumeElement(const std::array<Vector3, 8> &vertices, bool wire) {
//}

struct BoxInstanceData {
    Matrix4 transformationMatrix;
    Color3 color;
};

struct Boxes {

    Boxes(Trade::MeshData mdata = Primitives::cubeWireframe()) {
        boxShader = Shaders::FlatGL3D{
                    Shaders::FlatGL3D::Flag::VertexColor |
                    Shaders::FlatGL3D::Flag::InstancedTransformation };
        boxInstanceBuffer = GL::Buffer{};
        boxMesh = MeshTools::compile(mdata);
        boxMesh.addVertexBufferInstanced(boxInstanceBuffer, 1, 0,
            Shaders::FlatGL3D::TransformationMatrix{},
            Shaders::FlatGL3D::Color3{});
    }

    void DrawFrame(const Matrix4 &projectionMatrix) {
        boxInstanceBuffer.setData(data, GL::BufferUsage::DynamicDraw);
        boxMesh.setInstanceCount(data.size());
        boxShader.setTransformationProjectionMatrix(projectionMatrix)
            .draw(boxMesh);
    }

    void Append(const Matrix4& mat4, const Color3 &color) {
        arrayAppend(data, InPlaceInit, mat4, color);
    }

    void Clear() {
        arrayResize(data, 0);
    }

    GL::Mesh boxMesh{ NoCreate };
    GL::Buffer boxInstanceBuffer{ NoCreate };
    Shaders::FlatGL3D boxShader{ NoCreate };
    Containers::Array<BoxInstanceData> data;
};


#endif