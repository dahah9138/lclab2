#include "Torus.h"


namespace LC {

void Torus::Init() {

    using namespace Magnum;
    NX = 4;
    NY = 4;
    Float CX = 64.0f;
    Float CY = 64.0f;

    Float r = 1.0f;
    Float R = 2.0f;

    // Map [0,1]x[0,1] to [0,2Pi]x[0,2Pi]

    Float dTheta = 2.0f * M_PI / NX;
    Float dPhi = 2.0f * M_PI / NY;


    UnsignedInt iModx = NX - 1;
    UnsignedInt iMody = NY - 1;

    UnsignedInt numVerts = NX * NY;
    UnsignedInt numInds = 6 * NX * NY;

    data = Containers::Array<Vertex>{ NoInit, numVerts };
    std::vector<UnsignedInt> indVec;
    indVec.reserve(numInds);

    for (std::size_t i = 0; i < numVerts; ++i) {

        UnsignedInt ii = i / NY;
        UnsignedInt jj = i - ii * NY;
        Float theta = ii * dTheta;
        Float phi = jj * dPhi;

        // Replace with torus parametrization
        data[i].position = Vector3{ (r * cos(theta) + R) * cos(phi), (r * cos(theta) + R) * sin(phi), r * sin(theta) };

        Float r = data[i].position.length();

        data[i].color = Color3{ cos(CX * r) * cos(CX * r), sin(CY * r) * sin(CY * r), 0.0f };

        // i = ii * NY + jj

        UnsignedInt right = (ii < iModx) ? ii + 1 : 0;
        UnsignedInt down = (jj < iMody) ? jj + 1 : 0;

        UnsignedInt ind_down = ii * NY + down;
        UnsignedInt ind_right = right * NY + jj;
        UnsignedInt ind_down_right = right * NY + down;



        indVec.emplace_back(i);
        indVec.emplace_back(ind_down);
        indVec.emplace_back(ind_right);
        
        indVec.emplace_back(ind_down);
        indVec.emplace_back(ind_down_right);
        indVec.emplace_back(ind_right);
        

        
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

void Torus::Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) {

    sheetShader.setTransformationProjectionMatrix(projection * arcball->viewMatrix())
               .draw(sheetMesh);
}

}