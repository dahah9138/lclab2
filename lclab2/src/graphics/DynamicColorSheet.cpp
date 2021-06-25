#include "DynamicColorSheet.h"

namespace LC {

    void DynamicColorSheet::Init(PositionFunction pos) {
        using namespace Magnum;

        //Default position function
        if (pos == 0) {
            pos = [](Float x, Float y, Float z) { return Vector3{x, y, z}; };
        }


        UnsignedInt numVerts = NX * NY;
        UnsignedInt numInds = 6 * (NX - 1) * (NY - 1);

        data = Containers::Array<Vertex>{ NoInit, numVerts };
        std::vector<UnsignedInt> indVec;
        indVec.reserve(numInds);

        for (std::size_t i = 0; i < numVerts; ++i) {

            // Matlab indexing
            UnsignedInt jj = i / NX;
            UnsignedInt ii = i - jj * NX;
            Float x = (Float)ii / (NX - 1);
            Float y = (Float)jj / (NY - 1);

            data[i].position = pos(-0.5f + x, CY / CX * (-0.5f + y), 0.0f);

            

            Float r = data[i].position.length();

            data[i].color = Color4{ 1.0f, 1.0f, 1.0f, 0.5f };

            if (ii < NX - 1 && jj < NY - 1)
            {
                // i = ii * NY + jj

                indVec.emplace_back(i);
                indVec.emplace_back(i + 1);
                indVec.emplace_back(NX + i);

                indVec.emplace_back(i + 1);
                indVec.emplace_back(NX + i + 1);
                indVec.emplace_back(NX + i);
            }
        }

        MeshIndexType indType;
        UnsignedInt indStart, indEnd;

        std::tie(ind_data, indType, indStart, indEnd) = MeshTools::compressIndices(indVec);

        sheetShader = Shaders::VertexColorGL3D{};
        sheetBuffer = GL::Buffer{};
        sheetIndexBuffer = GL::Buffer{};

        sheetBuffer.setData(data, GL::BufferUsage::DynamicDraw);
        sheetIndexBuffer.setData(ind_data, GL::BufferUsage::StaticDraw);

        sheetMesh.setCount(numInds).addVertexBuffer(sheetBuffer, 0,
            Shaders::VertexColorGL3D::Position{},
            Shaders::VertexColorGL3D::Color4{})
            .setIndexBuffer(sheetIndexBuffer, 0, indType, indStart, indEnd);
    }

    void DynamicColorSheet::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {

        Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);

        sheetShader.setTransformationProjectionMatrix(projection * arcball->viewMatrix())
            .draw(sheetMesh);

        Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::FaceCulling);
    }


}