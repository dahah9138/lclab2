#include "DynamicColorSheet.h"

namespace LC {

    void DynamicColorSheet::Init(PositionFunction pos, Float offset) {

        //Default position function
        if (pos == 0) {
            pos = [](Float x, Float y, Float z) { return Vector3{ x, y, z }; };
        }


        UnsignedInt numVerts = NX * NY;
        UnsignedInt numInds = 6 * (NX - 1) * (NY - 1);

        vertices = Containers::Array<Vertex>{ NoInit, numVerts };
        indices = Containers::Array<UnsignedInt>{ NoInit, numInds };

        UnsignedInt indCount = 0;

        // Generate vertices and indices
        for (std::size_t i = 0; i < numVerts; ++i) {

            // Matlab indexing
            UnsignedInt jj = i / NX;
            UnsignedInt ii = i - jj * NX;
            Float x = (Float)ii / (NX - 1);
            Float y = (Float)jj / (NY - 1);

            vertices[i].position = pos(CX * (-0.5f + x), CY * (-0.5f + y), offset);

            Float r = vertices[i].position.length();

            vertices[i].color = Color4{ 1.0f, 1.0f, 1.0f, 0.5f };

            if (ii < NX - 1 && jj < NY - 1)
            {
                // i = ii * NY + jj
                indices[indCount++] = i;
                indices[indCount++] = i + 1;
                indices[indCount++] = NX + i;

                indices[indCount++] = i + 1;
                indices[indCount++] = NX + i + 1;
                indices[indCount++] = NX + i;
            }
        }


        // Set Attributes
        attributes = Containers::Array<Trade::MeshAttributeData>{ 2 };

        attributes[0] = Trade::MeshAttributeData{ Trade::MeshAttribute::Position,
            Containers::stridedArrayView(vertices, &vertices[0].position,
                Containers::arraySize(vertices), sizeof(Vertex)) };
        attributes[1] = Trade::MeshAttributeData{ Trade::MeshAttribute::Color,
            Containers::stridedArrayView(vertices, &vertices[0].color,
                Containers::arraySize(vertices), sizeof(Vertex)) };
    }

    Trade::MeshData DynamicColorSheet::Data() {
        return Trade::MeshData{ MeshPrimitive::Triangles,
            {}, indices, Trade::MeshIndexData{indices},
            {}, vertices, Trade::meshAttributeDataNonOwningArray(attributes) };
    }

    GL::Mesh DynamicColorSheet::Mesh() {
        Trade::MeshData data = Data();
        GL::Mesh mesh{ data.primitive() };

        vertexBuffer.setData(MeshTools::interleave(data.positions3DAsArray(), data.colorsAsArray()));
        mesh.addVertexBuffer(vertexBuffer, 0, Shaders::VertexColorGL3D::Position{}, Shaders::VertexColorGL3D::Color4{});
        // Set up index buffer
        if (data.isIndexed()) {

            std::pair<Containers::Array<char>, MeshIndexType> compressed = MeshTools::compressIndices(data.indicesAsArray());
            GL::Buffer indices;
            indices.setData(compressed.first);
            mesh.setIndexBuffer(std::move(indices), 0, compressed.second).setCount(data.indexCount());
        }

        return mesh;
    }

    void DynamicColorSheet::Set(UnsignedInt NX_, UnsignedInt NY_, Float CX_, Float CY_) {
        NX = NX_;
        NY = NY_;
        CX = CX_;
        CY = CY_;
    }


}