#include "Sheet.h"

namespace LC {

void Sheet::Init() {
    using namespace Magnum;
    NX = 400;
    NY = 400;
    Float CX = 64.0f;
    Float CY = 64.0f;

    UnsignedInt numVerts = NX * NY;
    UnsignedInt numInds = 6 * (NX - 1) * (NY - 1);
    UnsignedInt indCount = 0;

    data = Containers::Array<Vertex>{ NoInit, numVerts };
    std::vector<UnsignedInt> indVec;
    indVec.reserve(numInds);

    for (std::size_t i = 0; i < numVerts; ++i) {

        UnsignedInt ii = i / NY;
        UnsignedInt jj = i - ii * NY;
        Float x = (Float)ii / (NX - 1);
        Float y = (Float)jj / (NY - 1);

        data[i].position = Vector3(-0.5f + x, CY / CX * (-0.5f + y), 0.0f);

        Float r = data[i].position.length();

        data[i].color = Color3{ cos(CX * r) * cos(CX * r), sin(CY * r) * sin(CY * r), 0.0f };

        if (ii < NX - 1 && jj < NY - 1)
        {
            // i = ii * NY + jj

            indVec.emplace_back(i);
            indVec.emplace_back(i + 1);
            indVec.emplace_back(NY + i);

            indVec.emplace_back(i + 1);
            indVec.emplace_back(NY + i + 1);
            indVec.emplace_back(NY + i);
        }
    }

    MeshIndexType indType;
    UnsignedInt indStart, indEnd;

    std::tie(ind_data, indType, indStart, indEnd) = MeshTools::compressIndices(indVec);

    sheetShader = Shaders::VertexColorGL3D{};
    sheetBuffer = GL::Buffer{};
    sheetIndexBuffer = GL::Buffer{};

    sheetBuffer.setData(data, GL::BufferUsage::StaticDraw);
    sheetIndexBuffer.setData(ind_data, GL::BufferUsage::StaticDraw);

    sheetMesh.setCount(numInds).addVertexBuffer(sheetBuffer, 0,
        Shaders::VertexColorGL3D::Position{},
        Shaders::VertexColorGL3D::Color3{})
        .setIndexBuffer(sheetIndexBuffer, 0, indType, indStart, indEnd);
}

void Sheet::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {
    
    Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);

    sheetShader.setTransformationProjectionMatrix(projection * arcball->viewMatrix())
               .draw(sheetMesh);

    Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::FaceCulling);
}

}