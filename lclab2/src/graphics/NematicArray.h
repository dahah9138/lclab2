#ifndef NEMATICARRAY_H
#define NEMATICARRAY_H

#include "core.h"

#include <Magnum/ImageView.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Primitives/Capsule.h>
#include <Magnum/Primitives/Cone.h>
#include <Magnum/Primitives/Cylinder.h>
#include <Magnum/Primitives/Circle.h>
#include <Magnum/Trade/MeshData.h>
#include <Corrade/Containers/Optional.h>

#include "ArcBall.h"


namespace LC {

struct NematicArray {
	
    struct PolyInstanceData {
        Magnum::Matrix4 transformationMatrix;
        Magnum::Matrix3x3 normalMatrix;
        Magnum::Color3 color;
    };

    // Example of how multiple meshes can be compiled together
    // (Note texture coordinates would be necessary for the cylinder caps
    // to have different colors)
    struct Cylinder {
        static Magnum::GL::Mesh Compile(int rings, int segments, float halfLength, bool caps = false) {
            using namespace Magnum;
            // Create the vertex format used to compile
            struct MeshVertexData {
                Vector3 position;
                Vector3 normal;
            };

            Trade::MeshData cylinderBody = Primitives::cylinderSolid(rings, segments, halfLength);
            Trade::MeshData topCap = Primitives::circle3DSolid(segments);
            Trade::MeshData botCap = Primitives::circle3DSolid(segments);
            
            // Data stored as position, normal


            // Modify top cap and bot cap positions
            MeshVertexData* topCapVertexData =  (MeshVertexData*)&topCap.mutableVertexData()[0];
            MeshVertexData* botCapVertexData = (MeshVertexData*)&botCap.mutableVertexData()[0];
            for (int i = 0; i < topCap.vertexCount(); i++) {
                // Translate z positions
                topCapVertexData[i].position[2] += halfLength;
                botCapVertexData[i].position[2] -= halfLength;
                // Go from xy plane to xz plane
                // (So swap y and z)
                std::swap(topCapVertexData[i].position[1], topCapVertexData[i].position[2]);
                std::swap(topCapVertexData[i].normal[1], topCapVertexData[i].normal[2]);
                // Flip the normal of the bottom cap
                //botCapVertexData[i].normal[2] *= -1.;
                std::swap(botCapVertexData[i].position[1], botCapVertexData[i].position[2]);
                std::swap(botCapVertexData[i].normal[1], botCapVertexData[i].normal[2]);
            }

            // Combine everything into one mesh
            MeshVertexData* cylinderVertexData = (MeshVertexData*)&cylinderBody.mutableVertexData()[0];

            unsigned int* cylinderIndexData = (unsigned int*)&cylinderBody.mutableIndexData()[0];

            unsigned int numVerts = cylinderBody.vertexCount() + 2 * topCap.vertexCount();
            unsigned int numCircleInds = segments * 3;
            unsigned int numCircleVerts = segments + 1;
            unsigned int numInds = cylinderBody.indexCount() + 2 * numCircleInds;

            Containers::Array<Trade::MeshAttributeData> attributes;
            // Set Attributes
            attributes = Containers::Array<Trade::MeshAttributeData>{ 2 };
            Containers::Array<MeshVertexData> vertices = Containers::Array<MeshVertexData>{ NoInit, numVerts };;
            Containers::Array<UnsignedInt> indices = Containers::Array<UnsignedInt>{ NoInit, numInds };

            // Fill vertices
            {
                unsigned int counter = 0;
                for (unsigned int i = 0; i < cylinderBody.vertexCount(); i++) {
                    vertices[counter++] = cylinderVertexData[i];
                }
                for (unsigned int i = 0; i < topCap.vertexCount(); i++) {
                    vertices[counter++] = topCapVertexData[i];
                }
                for (unsigned int i = 0; i < botCap.vertexCount(); i++) {
                    vertices[counter++] = botCapVertexData[i];
                }
            }
            // Fill indices
            {
                unsigned int counter = 0;
                for (unsigned int i = 0; i < cylinderBody.indexCount(); i++) {
                    indices[counter++] = cylinderIndexData[i];
                }
                for (unsigned int i = 0; i < segments; i++) {
                    indices[counter++] = cylinderBody.indexCount();
                    unsigned int ip1 = (i + 1) % numCircleVerts;
                    unsigned int ip2 = (i + 2) % numCircleVerts;
                    if (ip1 == 0) ip1 = 1;
                    if (ip2 == 0) ip2 = 1;
                    indices[counter++] = ip1 + cylinderBody.indexCount();
                    indices[counter++] = ip2 + cylinderBody.indexCount();
                }
                for (unsigned int i = 0; i < segments; i++) {
                    indices[counter++] = cylinderBody.indexCount() + numCircleInds;
                    unsigned int ip1 = (i + 1) % numCircleVerts;
                    unsigned int ip2 = (i + 2) % numCircleVerts;
                    if (ip1 == 0) ip1 = 1;
                    if (ip2 == 0) ip2 = 1;
                    indices[counter++] = ip1 + cylinderBody.indexCount() + numCircleInds;
                    indices[counter++] = ip2 + cylinderBody.indexCount() + numCircleInds;
                }
            }

            attributes[0] = Trade::MeshAttributeData{ Trade::MeshAttribute::Position,
                Containers::stridedArrayView(vertices, &vertices[0].position,
                    Containers::arraySize(vertices), sizeof(MeshVertexData)) };
            attributes[1] = Trade::MeshAttributeData{ Trade::MeshAttribute::Normal,
                Containers::stridedArrayView(vertices, &vertices[0].normal,
                    Containers::arraySize(vertices), sizeof(MeshVertexData)) };

            Trade::MeshData data = Trade::MeshData{ MeshPrimitive::Triangles,
            {}, indices, Trade::MeshIndexData{indices},
            {}, vertices, Trade::meshAttributeDataNonOwningArray(attributes) };


            GL::Mesh mesh{ data.primitive() };

            // This is normally held onto...
            GL::Buffer vertexBuffer;
            vertexBuffer.setData(MeshTools::interleave(data.positions3DAsArray(), data.normalsAsArray()));
            mesh.addVertexBuffer(vertexBuffer, 0, Shaders::PhongGL::Position{}, Shaders::PhongGL::Normal{});
            // Set up index buffer
            if (data.isIndexed()) {
                std::pair<Containers::Array<char>, MeshIndexType> compressed = MeshTools::compressIndices(data.indicesAsArray());
                GL::Buffer indices;
                indices.setData(compressed.first);
                mesh.setIndexBuffer(std::move(indices), 0, compressed.second).setCount(data.indexCount());
            }

            return mesh;
        }

    };

    enum DrawType { Ellipsoid = 0, Cylinder = 1, Cone = 2};

    void Init(void* positions, std::function<Magnum::Vector3(void*, std::size_t)> Access, std::size_t size, int subdivisions = 2);
    void Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4 &projection);
    void Draw(const Magnum::Matrix4& viewMatrix, const Magnum::Matrix4& projection);

    static std::map<DrawType, std::string> DrawMap() { return { {Ellipsoid, "Ellipsoid" }, { Cylinder, "Cylinder" }, { Cone, "Cone" } }; }

    Magnum::UnsignedInt numObjects;

    Magnum::Containers::Array<Magnum::Vector3> polyPositions;
    Magnum::Float scale = 0.125f, hLength = 2.0f;

    Magnum::GL::Mesh polyMesh{ Magnum::NoCreate };

    std::vector< Magnum::GL::Mesh > polyMeshes;

    Magnum::GL::Buffer polyInstanceBuffer{ Magnum::NoCreate };
    Magnum::Shaders::PhongGL polyShader{ Magnum::NoCreate };
    Magnum::Containers::Array<PolyInstanceData> polyInstanceData;

    Magnum::Float specular = 0.f;
    Magnum::Float diffuse = 0.1f;
    Magnum::Float ambient = 1.f;
    std::map<DrawType, std::string> map = DrawMap();
    DrawType selected_drawType = DrawType::Cylinder;
};
}

#endif