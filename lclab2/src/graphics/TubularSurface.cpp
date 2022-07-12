#include "TubularSurface.h"

namespace LC {

void TubularSurface::Init(const bool &closedTube) {
    
	using namespace Magnum;
	
    NX = 32;
    // Number of curve points
    NY = curve.size();

    const Float r = radius;
    const Float R = 1.0f;

    // Map f(u,v): [0,1]x[0,1] -> [0,2Pi)x[0,2Pi)

    Float dTheta = 2.0f * M_PI / NX;


    UnsignedInt iModx = NX - 1;
    UnsignedInt iMody = NY - 1;

    UnsignedInt numVerts = NX * NY;
    UnsignedInt numInds = 6 * NX * (NY - !closedTube);

    // Invalid curve
    if (!numVerts) return;

    // Tangent vector
    // Input: current point of curve
    // Output: Tangent vector
    auto Tangent = [&](int ti) {
        // PBCs
        int tm1 = ti - 1;
        int tp1 = ti + 1;

        tp1 = tp1 % NY;
        if (tm1 < 0) tm1 = NY - 1;

        Eigen::Vector3d result = curve[tp1] - curve[tm1];
        result.normalize();
        return result;
    };

    // Binormal vector
    // Input: current point of curve
    // Output: Binormal vector
    auto Binormal = [&](int ti) {
        // PBCs
        int tm1 = ti - 1;
        int tp1 = ti + 1;

        tp1 = tp1 % NY;
        if (tm1 < 0) tm1 = NY - 1;

        Eigen::Vector3d v1 = curve[ti] - curve[tm1];
        Eigen::Vector3d v2 = curve[tp1] - curve[ti];
        Eigen::Vector3d result = v1.cross(v2);
        result.normalize();
        return result;
    };

    // Tangent vector
    // Input: current point of curve
    // Output: Tangent vector
    auto Normal = [&](int ti) {
        return Binormal(ti).cross(Tangent(ti));
    };

    // Tubular surface point given
    // u - Circle parameter about curve
    // v - Point on the curve
    auto TubePoint = [&](float u, int t) {
        // Compute TNB frame
        Eigen::Vector3d tangent = Tangent(t);
        Eigen::Vector3d binormal = Binormal(t);
        Eigen::Vector3d normal = binormal.cross(tangent);

        Eigen::Vector3d rho = r * cos(u) * normal + r * sin(u) * binormal;
        Eigen::Vector3d point = curve[t] + rho;

        // Make sure that the point returned is a Magnum vector point
        Vector3 mpoint, surface_normal;

        rho.normalize();

        for (int d = 0; d < 3; d++) {
            mpoint[d] = point(d);
            surface_normal[d] = rho(d);
        }

        std::array<Vector3, 2> pnorm = { mpoint, surface_normal };

        return pnorm;
    };

    data = Containers::Array<Vertex>{ NoInit, numVerts };
    indices.reserve(numInds);

    for (std::size_t i = 0; i < numVerts; ++i) {


        UnsignedInt t = i / NX;
        // Curve parametrization
        UnsignedInt jj = i - t * NX; // v
        Float u = jj * dTheta;

        auto pnorm_tuple = TubePoint(u, t);

        data[i].position = pnorm_tuple[0];
        data[i].normal = pnorm_tuple[1];

        // Color wheel
        Deg d(u * 180.0f / M_PI);
        data[i].color = Color4::fromHsv({ d, 1.0f, 1.0f }, alpha);
		

        UnsignedInt down = (t + 1) % NY;
        UnsignedInt right = (jj + 1) % NX;

        UnsignedInt ind_down = down * NX + jj;
        UnsignedInt ind_right = t * NX + right;
        UnsignedInt ind_down_right = down * NX + right;

        // Either tube is closed (PBCs) or t is one row from bottom ring
        if (closedTube || t < NY - 1) {

            indices.emplace_back(i);
            indices.emplace_back(ind_right);
            indices.emplace_back(ind_down);
        
            indices.emplace_back(ind_down);
            indices.emplace_back(ind_right);
            indices.emplace_back(ind_down_right);

        }
    }

    // Tube smoothing

    // Step 1. Untwist the tube by fixing ill-defined normals
    bool notSmooth = true;

    Eigen::Vector3d tangent, binormal, normal;
    tangent = Tangent(0);
    binormal = Binormal(0);
    normal = binormal.cross(tangent);

    for (int j = 0; j < NY - 1; j++) { // jth ring

        bool foundBadPoint = false;

        // Get Frenet frames for (j+1)th ring
        Eigen::Vector3d tangent2 = Tangent(j + 1);
        Eigen::Vector3d binormal2 = Binormal(j + 1);
        Eigen::Vector3d normal2 = binormal2.cross(tangent2);

        float dnormal = normal.dot(normal2);
        float dbinormal = binormal.dot(binormal2);

        for (int i = 0; i < NX; i++) {

            float u = i * dTheta;
            Eigen::Vector3d rho, pt;

            // If only normal flipped
            if (dnormal < 0.f && dbinormal > 0.f) {
                normal = -normal2;
                rho = r * cos(u) * -normal2 + r * sin(u) * binormal2;
                foundBadPoint = true;
            }
            // Only binormal flipped
            else if (dnormal > 0.f && dbinormal < 0.f) {
                binormal = -binormal2;
                rho = r * cos(u) * normal2 + r * sin(u) * -binormal2;
                foundBadPoint = true;
            }
            // Both flipped
            else if (dnormal < 0.f && dbinormal < 0.f) {
                normal = -normal2;
                binormal = -binormal2;
                rho = r * cos(u) * -normal2 + r * sin(u) * -binormal2;
                foundBadPoint = true;
            }
            else { // These two frenet frames are fine
                normal = normal2;
                binormal = binormal2;
                break;
            }

            pt = curve[j + 1] + rho;
            data[i + (j + 1) * NX].position = Vector3(pt[0], pt[1], pt[2]);
            rho.normalize();
            data[i + (j + 1) * NX].normal = Vector3(rho[0], rho[1], rho[2]);

        }

    }

    for (int steps = 0; steps < relaxIterations; steps++)
        // Step 2. Rotate each ring to an equilibrium position
        for (int j = 0; j < NY - !closedTube; j++) { // jth ring

            int jp1 = (j + 1) % NY;

            // Compute quaternion for jth ring
            Eigen::Vector3d tangent2 = Tangent(j + 1);
            Vector3 rho = data[0 + j * NX].position - Vector3(curve[j].x(), curve[j].y(), curve[j].z());
            Vector3 rho2 = data[0 + jp1 * NX].position - Vector3(curve[jp1].x(), curve[jp1].y(), curve[jp1].z());
            float alpha = 0.05f * M_PI;
            float sign = dot(cross(Vector3(tangent2.x(), tangent2.y(), tangent2.z()), rho2), rho - rho2);
            if (sign != 0.)
                sign /= abs(sign);
            float force = sign * alpha / (r * r) * cross(rho, rho2).length();

            tangent2 = sin(force / 2.) * tangent2;

            Quaternion quat(
                Vector3(tangent2.x(),
                tangent2.y(),
                tangent2.z()),
                cos(force / 2.)
            );

            auto rotMatrix = quat.toMatrix();

            // Update points for (j + 1)th ring
            for (int i = 0; i < NX; i++) {
                rho2 = data[i + jp1 * NX].position - Vector3(curve[jp1].x(), curve[jp1].y(), curve[jp1].z());
                rho2 = rotMatrix * rho2;
                data[i + jp1 * NX].position = rho2 + Vector3(curve[jp1].x(), curve[jp1].y(), curve[jp1].z());
                data[i + jp1 * NX].normal = rho2 / rho2.length();
            }
        }

    if (compile) {

        // Compress indices

        MeshIndexType indType;
        UnsignedInt indStart, indEnd;

        std::tie(ind_data, indType, indStart, indEnd) = MeshTools::compressIndices(indices);

        sheetShader = Shaders::PhongGL{ Shaders::PhongGL::Flag::VertexColor };
        sheetBuffer = GL::Buffer{};
        sheetIndexBuffer = GL::Buffer{};

        sheetBuffer.setData(data, GL::BufferUsage::StaticDraw);
        sheetIndexBuffer.setData(ind_data, GL::BufferUsage::StaticDraw);

        sheetMesh.setCount(numInds).addVertexBuffer(sheetBuffer, 0,
            Shaders::PhongGL::Position{},
            Shaders::PhongGL::Normal{},
		    Shaders::PhongGL::Color4{})
            .setIndexBuffer(sheetIndexBuffer, 0, indType, indStart, indEnd);
    }
}

}