#include "Surface.h"

namespace LC {

    void Surface::Init(Vertex *verts, unsigned int nVerts, unsigned int *inds, unsigned int nIndices) {

        vertices = Containers::Array<Vertex>{ NoInit, nVerts };
        indices = Containers::Array<UnsignedInt>{ NoInit, nIndices };

        // Generate vertices and indices
        for (std::size_t i = 0; i < nVerts; ++i) {
            vertices[i] = verts[i];
        }
        for (std::size_t i = 0; i < nIndices; ++i) {
            indices[i] = inds[i];
        }

        // Set Attributes
        attributes = Containers::Array<Trade::MeshAttributeData>{ 3 };

        attributes[0] = Trade::MeshAttributeData{ Trade::MeshAttribute::Position,
            Containers::stridedArrayView(vertices, &vertices[0].position,
                Containers::arraySize(vertices), sizeof(Vertex)) };
        attributes[1] = Trade::MeshAttributeData{ Trade::MeshAttribute::Normal,
            Containers::stridedArrayView(vertices, &vertices[0].normal,
                Containers::arraySize(vertices), sizeof(Vertex)) };
        attributes[2] = Trade::MeshAttributeData{ Trade::MeshAttribute::Color,
            Containers::stridedArrayView(vertices, &vertices[0].color,
                Containers::arraySize(vertices), sizeof(Vertex)) };
    }

    Trade::MeshData Surface::Data() {
        return Trade::MeshData{ MeshPrimitive::Triangles,
            {}, indices, Trade::MeshIndexData{indices},
            {}, vertices, Trade::meshAttributeDataNonOwningArray(attributes) };
    }

    GL::Mesh Surface::Mesh() {
        Trade::MeshData data = Data();
        GL::Mesh mesh{ data.primitive() };

        vertexBuffer.setData(MeshTools::interleave(data.positions3DAsArray(), data.normalsAsArray(), data.colorsAsArray()));
        mesh.addVertexBuffer(vertexBuffer, 0, Shaders::PhongGL::Position{}, Shaders::PhongGL::Normal{}, Shaders::PhongGL::Color4{});
        // Set up index buffer
        if (data.isIndexed()) {

            std::pair<Containers::Array<char>, MeshIndexType> compressed = MeshTools::compressIndices(data.indicesAsArray());
            GL::Buffer indices;
            indices.setData(compressed.first);
            mesh.setIndexBuffer(std::move(indices), 0, compressed.second).setCount(data.indexCount());
        }

        return mesh;
    }


}