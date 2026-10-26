#include "MeshProcessor.h"

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <vtkNew.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkFlyingEdges3D.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellArrayIterator.h>

namespace LC {
    namespace Math {
        struct MeshProcessor::Pipeline {
            vtkNew<vtkFloatArray> scalars;
            vtkNew<vtkImageData> volume;
            vtkNew<vtkFlyingEdges3D> contour;
        };

        MeshProcessor::MeshProcessor()
            : vtk(std::make_unique<Pipeline>()) {
            vtk->volume->SetSpacing(1.0, 1.0, 1.0);
            vtk->volume->SetOrigin(0.0, 0.0, 0.0);
            vtk->volume->GetPointData()->SetScalars(vtk->scalars);

            vtk->contour->SetInputData(vtk->volume);
            vtk->contour->SetValue(0, iso);
            vtk->contour->ComputeNormalsOn();
            vtk->contour->ComputeGradientsOff();
            vtk->contour->ComputeScalarsOff();
            vtk->contour->InterpolateAttributesOff();
        }

        MeshProcessor::MeshProcessor(int Nx, int Ny, int Nz, float isoValue)
            : MeshProcessor() {
            SetResolution(Nx, Ny, Nz);
            SetIsoValue(isoValue);
        }

        MeshProcessor::~MeshProcessor() = default;
        MeshProcessor::MeshProcessor(MeshProcessor&&) noexcept = default;
        MeshProcessor& MeshProcessor::operator=(MeshProcessor&&) noexcept = default;

        void MeshProcessor::SetResolution(int Nx, int Ny, int Nz) {
            if (nx == Nx && ny == Ny && nz == Nz) return;
            nx = Nx;
            ny = Ny;
            nz = Nz;
            vtk->volume->SetDimensions(nx, ny, nz);
        }

        void MeshProcessor::SetSpacing(float dx, float dy, float dz) {
            vtk->volume->SetSpacing(dx, dy, dz);
        }

        void MeshProcessor::SetOrigin(float ox, float oy, float oz) {
            vtk->volume->SetOrigin(ox, oy, oz);
        }

        void MeshProcessor::SetIsoValue(float isoValue) {
            iso = isoValue;
            vtk->contour->SetValue(0, iso);
        }

        MeshProcessor::ExtractedMesh MeshProcessor::ComputeIsosurface(const float* voxelData) {
            vtkIdType totalSize = static_cast<vtkIdType>(nx) * ny * nz;
            if (totalSize <= 0 || !voxelData) return {};

            vtk->scalars->SetArray(const_cast<float*>(voxelData), totalSize, 1);
            vtk->scalars->Modified();
            vtk->volume->Modified();

            vtk->contour->Update();
            return ExtractPolyData(vtk->contour->GetOutput(), true);
        }

        MeshProcessor::ExtractedMesh MeshProcessor::ComputeIsosurface(const float* voxelData, float isoValue) {
            SetIsoValue(isoValue);
            return ComputeIsosurface(voxelData);
        }

        void MeshProcessor::BuildOrthonormalBasis(
            const std::array<float, 3>& n0,
            std::array<float, 3>& b1,
            std::array<float, 3>& b2
        ) {
            float norm = std::sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
            std::array<float, 3> u = { n0[0] / norm, n0[1] / norm, n0[2] / norm };

            float sign = std::copysign(1.0f, u[2]);
            const float a = -1.0f / (sign + u[2]);
            const float b = u[0] * u[1] * a;

            b1 = { 1.0f + sign * u[0] * u[0] * a, sign * b, -sign * u[0] };
            b2 = { b, sign + u[1] * u[1] * a, -u[1] };
        }

        float MeshProcessor::SampleAlignmentTrilinear(
            const float* directorField,
            const std::array<float, 3>& n0,
            float x, float y, float z
        ) const {
            double origin[3], spacing[3];
            vtk->volume->GetOrigin(origin);
            vtk->volume->GetSpacing(spacing);

            float gx = std::clamp((x - static_cast<float>(origin[0])) / static_cast<float>(spacing[0]), 0.0f, static_cast<float>(nx - 1));
            float gy = std::clamp((y - static_cast<float>(origin[1])) / static_cast<float>(spacing[1]), 0.0f, static_cast<float>(ny - 1));
            float gz = std::clamp((z - static_cast<float>(origin[2])) / static_cast<float>(spacing[2]), 0.0f, static_cast<float>(nz - 1));

            int i0 = std::clamp(static_cast<int>(gx), 0, nx - 2);
            int j0 = std::clamp(static_cast<int>(gy), 0, ny - 2);
            int k0 = std::clamp(static_cast<int>(gz), 0, nz - 2);

            float tx = gx - i0;
            float ty = gy - j0;
            float tz = gz - k0;

            size_t strideY = static_cast<size_t>(nx);
            size_t strideZ = strideY * ny;
            size_t vol = strideZ * nz;
            size_t base = static_cast<size_t>(i0) + strideY * j0 + strideZ * k0;

            const float* nxPtr = directorField;
            const float* nyPtr = directorField + vol;
            const float* nzPtr = directorField + 2 * vol;

            auto dotAtOffset = [&](size_t idx) {
                return nxPtr[idx] * n0[0] + nyPtr[idx] * n0[1] + nzPtr[idx] * n0[2];
                };

            float c00 = dotAtOffset(base) * (1.0f - tx) + dotAtOffset(base + 1) * tx;
            float c10 = dotAtOffset(base + strideY) * (1.0f - tx) + dotAtOffset(base + strideY + 1) * tx;
            float c01 = dotAtOffset(base + strideZ) * (1.0f - tx) + dotAtOffset(base + strideZ + 1) * tx;
            float c11 = dotAtOffset(base + strideZ + strideY) * (1.0f - tx) + dotAtOffset(base + strideZ + strideY + 1) * tx;

            float c0 = c00 * (1.0f - ty) + c10 * ty;
            float c1 = c01 * (1.0f - ty) + c11 * ty;

            return c0 * (1.0f - tz) + c1 * tz;
        }

        MeshProcessor::ExtractedLines MeshProcessor::ComputeBasisPreimageLines(
            const float* directorField,
            const std::array<float, 3>& n0,
            float alignmentThreshold
        ) {
            ExtractedLines result;
            size_t vol = static_cast<size_t>(nx) * ny * nz;
            if (vol == 0 || !directorField || nx < 2 || ny < 2 || nz < 2) return result;

            float norm = std::sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
            if (norm < 1e-8f) return result;
            std::array<float, 3> u = { n0[0] / norm, n0[1] / norm, n0[2] / norm };

            std::array<float, 3> b1, b2;
            BuildOrthonormalBasis(u, b1, b2);

            double originD[3], spacingD[3];
            vtk->volume->GetOrigin(originD);
            vtk->volume->GetSpacing(spacingD);
            const float ox = static_cast<float>(originD[0]);
            const float oy = static_cast<float>(originD[1]);
            const float oz = static_cast<float>(originD[2]);
            const float dx = static_cast<float>(spacingD[0]);
            const float dy = static_cast<float>(spacingD[1]);
            const float dz = static_cast<float>(spacingD[2]);

            const float* nxPtr = directorField;
            const float* nyPtr = directorField + vol;
            const float* nzPtr = directorField + 2 * vol;
            const size_t strideY = static_cast<size_t>(nx);
            const size_t strideZ = strideY * ny;

            auto intersectTri = [&](
                const float pA[3], float f1A, float f2A, float alA,
                const float pB[3], float f1B, float f2B, float alB,
                const float pC[3], float f1C, float f2C, float alC,
                float outPt[3]
                ) -> bool {
                    if ((f1A > 0.0f && f1B > 0.0f && f1C > 0.0f) || (f1A < 0.0f && f1B < 0.0f && f1C < 0.0f) ||
                        (f2A > 0.0f && f2B > 0.0f && f2C > 0.0f) || (f2A < 0.0f && f2B < 0.0f && f2C < 0.0f)) {
                        return false;
                    }
                    if (alA < 0.0f && alB < 0.0f && alC < 0.0f) {
                        return false;
                    }

                    float d1B = f1B - f1A, d1C = f1C - f1A;
                    float d2B = f2B - f2A, d2C = f2C - f2A;
                    float det = d1B * d2C - d1C * d2B;
                    if (std::abs(det) < 1e-12f) return false;

                    float invDet = 1.0f / det;
                    float uB = (-f1A * d2C + f2A * d1C) * invDet;
                    float vC = (-d1B * f2A + d2B * f1A) * invDet;
                    float wA = 1.0f - uB - vC;

                    const float eps = -1e-5f;
                    if (uB < eps || vC < eps || wA < eps) return false;

                    float align = wA * alA + uB * alB + vC * alC;
                    if (align < alignmentThreshold) return false;

                    outPt[0] = wA * pA[0] + uB * pB[0] + vC * pC[0];
                    outPt[1] = wA * pA[1] + uB * pB[1] + vC * pC[1];
                    outPt[2] = wA * pA[2] + uB * pB[2] + vC * pC[2];
                    return true;
                };

            static const int faceQuads[6][4] = {
                { 0, 2, 6, 4 },
                { 1, 3, 7, 5 },
                { 0, 1, 5, 4 },
                { 2, 3, 7, 6 },
                { 0, 1, 3, 2 },
                { 4, 5, 7, 6 }
            };

            for (int k = 0; k < nz - 1; ++k) {
                for (int j = 0; j < ny - 1; ++j) {
                    size_t baseRow = strideY * j + strideZ * k;
                    for (int i = 0; i < nx - 1; ++i) {
                        size_t c000 = baseRow + i;
                        size_t cIdx[8] = {
                            c000,
                            c000 + 1,
                            c000 + strideY,
                            c000 + strideY + 1,
                            c000 + strideZ,
                            c000 + strideZ + 1,
                            c000 + strideZ + strideY,
                            c000 + strideZ + strideY + 1
                        };

                        float f1[8], f2[8], al[8];
                        int maskPos1 = 0, maskNeg1 = 0, maskPos2 = 0, maskNeg2 = 0, maskAlign = 0;

                        for (int c = 0; c < 8; ++c) {
                            size_t idx = cIdx[c];
                            float vx = nxPtr[idx], vy = nyPtr[idx], vz = nzPtr[idx];
                            al[c] = vx * u[0] + vy * u[1] + vz * u[2];
                            if (al[c] >= alignmentThreshold) maskAlign = 1;

                            f1[c] = vx * b1[0] + vy * b1[1] + vz * b1[2];
                            f2[c] = vx * b2[0] + vy * b2[1] + vz * b2[2];

                            if (f1[c] >= 0.0f) maskPos1 = 1; else maskNeg1 = 1;
                            if (f2[c] >= 0.0f) maskPos2 = 1; else maskNeg2 = 1;
                        }

                        if (!maskAlign || !maskPos1 || !maskNeg1 || !maskPos2 || !maskNeg2) {
                            continue;
                        }

                        float x0 = ox + i * dx, x1 = x0 + dx;
                        float y0 = oy + j * dy, y1 = y0 + dy;
                        float z0 = oz + k * dz, z1 = z0 + dz;

                        const float pos[8][3] = {
                            { x0, y0, z0 },{ x1, y0, z0 },{ x0, y1, z0 },{ x1, y1, z0 },
                            { x0, y0, z1 },{ x1, y0, z1 },{ x0, y1, z1 },{ x1, y1, z1 }
                        };

                        float hitPts[6][3];
                        int hitCount = 0;

                        for (int f = 0; f < 6; ++f) {
                            int v0 = faceQuads[f][0], v1 = faceQuads[f][1], v2 = faceQuads[f][2], v3 = faceQuads[f][3];
                            float pt[3];
                            if (intersectTri(pos[v0], f1[v0], f2[v0], al[v0],
                                pos[v1], f1[v1], f2[v1], al[v1],
                                pos[v2], f1[v2], f2[v2], al[v2], pt) ||
                                intersectTri(pos[v0], f1[v0], f2[v0], al[v0],
                                    pos[v2], f1[v2], f2[v2], al[v2],
                                    pos[v3], f1[v3], f2[v3], al[v3], pt)) {

                                bool duplicate = false;
                                for (int h = 0; h < hitCount; ++h) {
                                    float ddx = pt[0] - hitPts[h][0];
                                    float ddy = pt[1] - hitPts[h][1];
                                    float ddz = pt[2] - hitPts[h][2];
                                    if (ddx * ddx + ddy * ddy + ddz * ddz < 1e-10f) {
                                        duplicate = true;
                                        break;
                                    }
                                }
                                if (!duplicate && hitCount < 6) {
                                    hitPts[hitCount][0] = pt[0];
                                    hitPts[hitCount][1] = pt[1];
                                    hitPts[hitCount][2] = pt[2];
                                    ++hitCount;
                                }
                            }
                        }

                        if (hitCount >= 2) {
                            unsigned int baseIdx = static_cast<unsigned int>(result.vertices.size() / 3);
                            result.vertices.push_back(hitPts[0][0]);
                            result.vertices.push_back(hitPts[0][1]);
                            result.vertices.push_back(hitPts[0][2]);
                            result.vertices.push_back(hitPts[1][0]);
                            result.vertices.push_back(hitPts[1][1]);
                            result.vertices.push_back(hitPts[1][2]);
                            result.lineIndices.push_back(baseIdx);
                            result.lineIndices.push_back(baseIdx + 1);
                        }
                    }
                }
            }

            return result;
        }

        MeshProcessor::ExtractedMesh MeshProcessor::BuildTubeFromExtractedLines(const ExtractedLines& lines, float tubeRadius, int tubeSides) {
            if (lines.empty() || tubeSides < 3) return {};

            double spacing[3];
            vtk->volume->GetSpacing(spacing);
            float minStep = static_cast<float>(std::min({ spacing[0], spacing[1], spacing[2] }));

            using Pt3 = std::array<float, 3>;
            auto distSq = [](const Pt3& a, const Pt3& b) {
                float dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
                return dx * dx + dy * dy + dz * dz;
                };

            // Spatial hashing to weld vertices within 0.35 * minStep
            float weldTol = 0.35f * minStep;
            float invCell = (weldTol > 1e-8f) ? (1.0f / weldTol) : 1e6f;

            struct GridKey {
                int32_t x, y, z;
                bool operator==(const GridKey& o) const { return x == o.x && y == o.y && z == o.z; }
            };
            struct GridKeyHash {
                size_t operator()(const GridKey& k) const {
                    return (static_cast<size_t>(k.x) * 73856093u) ^
                        (static_cast<size_t>(k.y) * 19349663u) ^
                        (static_cast<size_t>(k.z) * 83492791u);
                }
            };

            std::unordered_map<GridKey, std::vector<uint32_t>, GridKeyHash> cellMap;
            std::vector<Pt3> weldedPts;
            size_t nRawVerts = lines.numVertices();
            weldedPts.reserve(nRawVerts / 2 + 1);
            std::vector<uint32_t> remap(nRawVerts);

            float weldTolSq = weldTol * weldTol;
            for (size_t i = 0; i < nRawVerts; ++i) {
                Pt3 p = { lines.vertices[3 * i], lines.vertices[3 * i + 1], lines.vertices[3 * i + 2] };
                int32_t gx = static_cast<int32_t>(std::floor(p[0] * invCell));
                int32_t gy = static_cast<int32_t>(std::floor(p[1] * invCell));
                int32_t gz = static_cast<int32_t>(std::floor(p[2] * invCell));

                uint32_t foundId = UINT32_MAX;
                float bestD2 = weldTolSq;

                for (int dz = -1; dz <= 1; ++dz) {
                    for (int dy = -1; dy <= 1; ++dy) {
                        for (int dx = -1; dx <= 1; ++dx) {
                            auto it = cellMap.find({ gx + dx, gy + dy, gz + dz });
                            if (it == cellMap.end()) continue;
                            for (uint32_t cand : it->second) {
                                float d2 = distSq(p, weldedPts[cand]);
                                if (d2 <= bestD2) {
                                    bestD2 = d2;
                                    foundId = cand;
                                }
                            }
                        }
                    }
                }

                if (foundId == UINT32_MAX) {
                    foundId = static_cast<uint32_t>(weldedPts.size());
                    weldedPts.push_back(p);
                    cellMap[{ gx, gy, gz }].push_back(foundId);
                }
                remap[i] = foundId;
            }

            // Build adjacency list and strip contiguous segments into polyline chains
            struct Edge { uint32_t to; size_t edgeIdx; };
            std::vector<std::vector<Edge>> adj(weldedPts.size());
            size_t numEdges = 0;
            for (size_t i = 0; i + 1 < lines.lineIndices.size(); i += 2) {
                uint32_t u = remap[lines.lineIndices[i]];
                uint32_t v = remap[lines.lineIndices[i + 1]];
                if (u != v) {
                    adj[u].push_back({ v, numEdges });
                    adj[v].push_back({ u, numEdges });
                    ++numEdges;
                }
            }
            if (numEdges == 0) return {};

            std::vector<uint8_t> visitedEdge(numEdges, 0);
            std::vector<std::vector<Pt3>> chains;

            auto traceFrom = [&](uint32_t startNode) {
                for (const auto& e : adj[startNode]) {
                    if (visitedEdge[e.edgeIdx]) continue;
                    visitedEdge[e.edgeIdx] = 1;

                    std::vector<Pt3> chain;
                    chain.push_back(weldedPts[startNode]);
                    uint32_t curr = e.to;
                    chain.push_back(weldedPts[curr]);

                    while (true) {
                        uint32_t nextNode = UINT32_MAX;
                        for (const auto& ne : adj[curr]) {
                            if (!visitedEdge[ne.edgeIdx]) {
                                visitedEdge[ne.edgeIdx] = 1;
                                nextNode = ne.to;
                                break;
                            }
                        }
                        if (nextNode == UINT32_MAX) break;
                        if (distSq(chain.back(), weldedPts[nextNode]) > 1e-8f) {
                            chain.push_back(weldedPts[nextNode]);
                        }
                        curr = nextNode;
                    }
                    if (chain.size() >= 2) {
                        chains.push_back(std::move(chain));
                    }
                }
                };

            // Trace from non-manifold or open endpoints first, then closed loops
            for (uint32_t i = 0; i < weldedPts.size(); ++i) {
                if (adj[i].size() != 2 && !adj[i].empty()) traceFrom(i);
            }
            for (uint32_t i = 0; i < weldedPts.size(); ++i) {
                if (adj[i].size() == 2) traceFrom(i);
            }

            if (chains.size() > 1) {
                chains.erase(
                    std::remove_if(chains.begin(), chains.end(), [](const std::vector<Pt3>& c) {
                        return c.size() < 4;
                        }),
                    chains.end()
                );
            }

            const float maxBridgeDistSq = (4.0f * minStep) * (4.0f * minStep);
            bool merged = true;
            while (merged && chains.size() > 1) {
                merged = false;
                float bestDist = maxBridgeDistSq;
                size_t bestI = 0, bestJ = 0;
                int bestMode = -1;

                for (size_t i = 0; i < chains.size(); ++i) {
                    if (chains[i].size() > 4 && distSq(chains[i].front(), chains[i].back()) < minStep * minStep) {
                        continue;
                    }
                    for (size_t j = i + 1; j < chains.size(); ++j) {
                        if (chains[j].size() > 4 && distSq(chains[j].front(), chains[j].back()) < minStep * minStep) {
                            continue;
                        }

                        float d0 = distSq(chains[i].back(), chains[j].front());
                        if (d0 < bestDist) { bestDist = d0; bestI = i; bestJ = j; bestMode = 0; }

                        float d1 = distSq(chains[i].back(), chains[j].back());
                        if (d1 < bestDist) { bestDist = d1; bestI = i; bestJ = j; bestMode = 1; }

                        float d2 = distSq(chains[i].front(), chains[j].front());
                        if (d2 < bestDist) { bestDist = d2; bestI = i; bestJ = j; bestMode = 2; }

                        float d3 = distSq(chains[i].front(), chains[j].back());
                        if (d3 < bestDist) { bestDist = d3; bestI = i; bestJ = j; bestMode = 3; }
                    }
                }

                if (bestMode != -1) {
                    auto& cA = chains[bestI];
                    auto& cB = chains[bestJ];
                    if (bestMode == 1) {
                        std::reverse(cB.begin(), cB.end());
                    }
                    else if (bestMode == 2) {
                        std::reverse(cA.begin(), cA.end());
                    }
                    else if (bestMode == 3) {
                        std::reverse(cA.begin(), cA.end());
                        std::reverse(cB.begin(), cB.end());
                    }

                    size_t startIdx = (distSq(cA.back(), cB.front()) < 0.05f * minStep * minStep) ? 1 : 0;
                    cA.insert(cA.end(), cB.begin() + startIdx, cB.end());
                    chains.erase(chains.begin() + bestJ);
                    merged = true;
                }
            }

            // Precompute unit circle table for tube cross-section
            std::vector<float> cosTable(tubeSides), sinTable(tubeSides);
            const float twoPi = 6.283185307179586f;
            for (int s = 0; s < tubeSides; ++s) {
                float angle = twoPi * static_cast<float>(s) / static_cast<float>(tubeSides);
                cosTable[s] = std::cos(angle);
                sinTable[s] = std::sin(angle);
            }

            auto normalize3 = [](Pt3& v) {
                float len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                if (len > 1e-8f) { v[0] /= len; v[1] /= len; v[2] /= len; }
                };
            auto cross3 = [](const Pt3& a, const Pt3& b) -> Pt3 {
                return {
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0]
                };
                };
            auto dot3 = [](const Pt3& a, const Pt3& b) -> float {
                return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                };

            ExtractedMesh result;
            const int smoothIters = 25;
            const float lambda = 0.5f;
            const float mu = -0.53f;
            const float loopSnapDistSq = (4.0f * minStep) * (4.0f * minStep);

            for (auto& curve : chains) {
                if (curve.size() < 3) continue;

                bool isClosed = (distSq(curve.front(), curve.back()) < loopSnapDistSq);
                if (isClosed && distSq(curve.front(), curve.back()) < 0.25f * minStep * minStep) {
                    curve.pop_back();
                }

                size_t count = curve.size();
                if (count < 3) continue;

                std::vector<Pt3> subdivided;
                subdivided.reserve(count * 2);
                for (size_t i = 0; i < count; ++i) {
                    subdivided.push_back(curve[i]);
                    if (isClosed || i + 1 < count) {
                        const Pt3& nextPt = curve[(i + 1) % count];
                        subdivided.push_back({
                            0.5f * (curve[i][0] + nextPt[0]),
                            0.5f * (curve[i][1] + nextPt[1]),
                            0.5f * (curve[i][2] + nextPt[2])
                            });
                    }
                }
                curve.swap(subdivided);
                count = curve.size();

                std::vector<Pt3> temp = curve;
                auto laplacianStep = [&](float weight) {
                    for (size_t i = 0; i < count; ++i) {
                        if (!isClosed && (i == 0 || i == count - 1)) {
                            temp[i] = curve[i];
                            continue;
                        }
                        size_t prev = (i == 0) ? (count - 1) : (i - 1);
                        size_t next = (i + 1) % count;
                        for (int d = 0; d < 3; ++d) {
                            float lap = 0.5f * (curve[prev][d] + curve[next][d]) - curve[i][d];
                            temp[i][d] = curve[i][d] + weight * lap;
                        }
                    }
                    curve.swap(temp);
                    };

                for (int iter = 0; iter < smoothIters; ++iter) {
                    laplacianStep(lambda);
                    laplacianStep(mu);
                }

                // Compute vertex tangents
                std::vector<Pt3> tangents(count);
                for (size_t i = 0; i < count; ++i) {
                    size_t prev = isClosed ? ((i == 0) ? count - 1 : i - 1) : (i == 0 ? 0 : i - 1);
                    size_t next = isClosed ? ((i + 1) % count) : (i + 1 < count ? i + 1 : count - 1);
                    tangents[i] = {
                        curve[next][0] - curve[prev][0],
                        curve[next][1] - curve[prev][1],
                        curve[next][2] - curve[prev][2]
                    };
                    normalize3(tangents[i]);
                }

                // Propagate parallel transport (Bishop) frame along curve
                std::vector<Pt3> normals(count), binormals(count);
                BuildOrthonormalBasis(tangents[0], normals[0], binormals[0]);

                for (size_t i = 1; i < count; ++i) {
                    const Pt3& Nprev = normals[i - 1];
                    const Pt3& Ti = tangents[i];
                    float proj = dot3(Nprev, Ti);
                    Pt3 Ncurr = { Nprev[0] - proj * Ti[0], Nprev[1] - proj * Ti[1], Nprev[2] - proj * Ti[2] };
                    normalize3(Ncurr);
                    normals[i] = Ncurr;
                    binormals[i] = cross3(Ti, Ncurr);
                }

                // For closed loops, untwist holonomy seam angle uniformly across all rings
                if (isClosed) {
                    const Pt3& Nlast = normals[count - 1];
                    const Pt3& T0 = tangents[0];
                    float proj = dot3(Nlast, T0);
                    Pt3 Ntransported = { Nlast[0] - proj * T0[0], Nlast[1] - proj * T0[1], Nlast[2] - proj * T0[2] };
                    normalize3(Ntransported);

                    float cosHol = std::clamp(dot3(Ntransported, normals[0]), -1.0f, 1.0f);
                    float sinHol = std::clamp(dot3(Ntransported, binormals[0]), -1.0f, 1.0f);
                    float holonomy = std::atan2(sinHol, cosHol);

                    for (size_t i = 0; i < count; ++i) {
                        float phi = -holonomy * static_cast<float>(i) / static_cast<float>(count);
                        float c = std::cos(phi), s = std::sin(phi);
                        Pt3 Ni = normals[i], Bi = binormals[i];
                        normals[i] = { c * Ni[0] + s * Bi[0], c * Ni[1] + s * Bi[1], c * Ni[2] + s * Bi[2] };
                        binormals[i] = { -s * Ni[0] + c * Bi[0], -s * Ni[1] + c * Bi[1], -s * Ni[2] + c * Bi[2] };
                    }
                }

                // Emit tube vertices and normals directly
                unsigned int baseVert = static_cast<unsigned int>(result.vertices.size());
                result.vertices.resize(baseVert + count * tubeSides);

                for (size_t i = 0; i < count; ++i) {
                    const Pt3& C = curve[i];
                    const Pt3& N = normals[i];
                    const Pt3& B = binormals[i];
                    for (int s = 0; s < tubeSides; ++s) {
                        float ca = cosTable[s], sa = sinTable[s];
                        float nx = ca * N[0] + sa * B[0];
                        float ny = ca * N[1] + sa * B[1];
                        float nz = ca * N[2] + sa * B[2];

                        auto& v = result.vertices[baseVert + i * tubeSides + s];
                        v.position[0] = C[0] + tubeRadius * nx;
                        v.position[1] = C[1] + tubeRadius * ny;
                        v.position[2] = C[2] + tubeRadius * nz;
                        v.normal[0] = nx;
                        v.normal[1] = ny;
                        v.normal[2] = nz;
                    }
                }

                // Emit tube triangle indices
                size_t numSegs = isClosed ? count : (count - 1);
                result.indices.reserve(result.indices.size() + numSegs * tubeSides * 6);

                for (size_t i = 0; i < numSegs; ++i) {
                    size_t nextI = (i + 1) % count;
                    unsigned int ring0 = baseVert + static_cast<unsigned int>(i * tubeSides);
                    unsigned int ring1 = baseVert + static_cast<unsigned int>(nextI * tubeSides);

                    for (int s = 0; s < tubeSides; ++s) {
                        int nextS = (s + 1) % tubeSides;
                        unsigned int i00 = ring0 + s;
                        unsigned int i01 = ring0 + nextS;
                        unsigned int i10 = ring1 + s;
                        unsigned int i11 = ring1 + nextS;

                        result.indices.push_back(i00);
                        result.indices.push_back(i10);
                        result.indices.push_back(i11);

                        result.indices.push_back(i00);
                        result.indices.push_back(i11);
                        result.indices.push_back(i01);
                    }
                }

                if (!isClosed) {
                    unsigned int capStartBase = static_cast<unsigned int>(result.vertices.size());
                    result.vertices.resize(capStartBase + 2 * (tubeSides + 1));

                    const Pt3& C0 = curve[0];
                    const Pt3& N0 = normals[0];
                    const Pt3& B0 = binormals[0];
                    const Pt3& T0 = tangents[0];

                    auto& startCenter = result.vertices[capStartBase];
                    startCenter.position[0] = C0[0];
                    startCenter.position[1] = C0[1];
                    startCenter.position[2] = C0[2];
                    startCenter.normal[0] = -T0[0];
                    startCenter.normal[1] = -T0[1];
                    startCenter.normal[2] = -T0[2];

                    for (int s = 0; s < tubeSides; ++s) {
                        float ca = cosTable[s], sa = sinTable[s];
                        float rx = ca * N0[0] + sa * B0[0];
                        float ry = ca * N0[1] + sa * B0[1];
                        float rz = ca * N0[2] + sa * B0[2];

                        auto& v = result.vertices[capStartBase + 1 + s];
                        v.position[0] = C0[0] + tubeRadius * rx;
                        v.position[1] = C0[1] + tubeRadius * ry;
                        v.position[2] = C0[2] + tubeRadius * rz;
                        v.normal[0] = -T0[0];
                        v.normal[1] = -T0[1];
                        v.normal[2] = -T0[2];
                    }

                    unsigned int capEndBase = capStartBase + tubeSides + 1;
                    size_t last = count - 1;
                    const Pt3& CL = curve[last];
                    const Pt3& NL = normals[last];
                    const Pt3& BL = binormals[last];
                    const Pt3& TL = tangents[last];

                    auto& endCenter = result.vertices[capEndBase];
                    endCenter.position[0] = CL[0];
                    endCenter.position[1] = CL[1];
                    endCenter.position[2] = CL[2];
                    endCenter.normal[0] = TL[0];
                    endCenter.normal[1] = TL[1];
                    endCenter.normal[2] = TL[2];

                    for (int s = 0; s < tubeSides; ++s) {
                        float ca = cosTable[s], sa = sinTable[s];
                        float rx = ca * NL[0] + sa * BL[0];
                        float ry = ca * NL[1] + sa * BL[1];
                        float rz = ca * NL[2] + sa * BL[2];

                        auto& v = result.vertices[capEndBase + 1 + s];
                        v.position[0] = CL[0] + tubeRadius * rx;
                        v.position[1] = CL[1] + tubeRadius * ry;
                        v.position[2] = CL[2] + tubeRadius * rz;
                        v.normal[0] = TL[0];
                        v.normal[1] = TL[1];
                        v.normal[2] = TL[2];
                    }

                    result.indices.reserve(result.indices.size() + tubeSides * 6);
                    for (int s = 0; s < tubeSides; ++s) {
                        int nextS = (s + 1) % tubeSides;

                        // Start cap (normal points along -T0, clockwise winding looking down +T0)
                        result.indices.push_back(capStartBase);
                        result.indices.push_back(capStartBase + 1 + nextS);
                        result.indices.push_back(capStartBase + 1 + s);

                        // End cap (normal points along +TL, counter-clockwise winding looking down +TL)
                        result.indices.push_back(capEndBase);
                        result.indices.push_back(capEndBase + 1 + s);
                        result.indices.push_back(capEndBase + 1 + nextS);
                    }
                }
            }

            return result;
        }

        MeshProcessor::ExtractedMesh MeshProcessor::ComputeBasisPreimageTube(
            const float* directorField,
            const std::array<float, 3>& n0,
            float tubeRadius,
            int tubeSides,
            float alignmentThreshold
        ) {
            ExtractedLines lines = ComputeBasisPreimageLines(directorField, n0, alignmentThreshold);
            return BuildTubeFromExtractedLines(lines, tubeRadius, tubeSides);
        }

        MeshProcessor::ExtractedLines MeshProcessor::ComputeUmbilicLines(
            const float* directorField,
            float alignmentThreshold
        ) {
            ExtractedLines result;
            size_t vol = static_cast<size_t>(nx) * ny * nz;
            if (vol == 0 || !directorField || nx < 3 || ny < 3 || nz < 3) return result;

            double originD[3], spacingD[3];
            vtk->volume->GetOrigin(originD);
            vtk->volume->GetSpacing(spacingD);
            const float ox = static_cast<float>(originD[0]);
            const float oy = static_cast<float>(originD[1]);
            const float oz = static_cast<float>(originD[2]);
            const float dx = static_cast<float>(spacingD[0]);
            const float dy = static_cast<float>(spacingD[1]);
            const float dz = static_cast<float>(spacingD[2]);

            const float* nxPtr = directorField;
            const float* nyPtr = directorField + vol;
            const float* nzPtr = directorField + 2 * vol;
            const size_t strideY = static_cast<size_t>(nx);
            const size_t strideZ = strideY * ny;

            // Precompute central-difference spatial gradients dn_b / dx_a at all grid nodes
            // Layout: 9 floats per voxel [dnx/dx, dny/dx, dnz/dx, dnx/dy, dny/dy, dnz/dy, dnx/dz, dny/dz, dnz/dz]
            std::vector<float> gradField(vol * 9, 0.0f);
            const float inv2dx = 0.5f / dx, invdx = 1.0f / dx;
            const float inv2dy = 0.5f / dy, invdy = 1.0f / dy;
            const float inv2dz = 0.5f / dz, invdz = 1.0f / dz;

            for (int k = 0; k < nz; ++k) {
                for (int j = 0; j < ny; ++j) {
                    size_t row = strideY * j + strideZ * k;
                    for (int i = 0; i < nx; ++i) {
                        size_t idx = row + i;
                        float nC[3] = { nxPtr[idx], nyPtr[idx], nzPtr[idx] };
                        float* g = &gradField[idx * 9];

                        auto alignedDiff = [&](size_t idxPlus, size_t idxMinus, float scale, float out[3]) {
                            float nP[3] = { nxPtr[idxPlus], nyPtr[idxPlus], nzPtr[idxPlus] };
                            float nM[3] = { nxPtr[idxMinus], nyPtr[idxMinus], nzPtr[idxMinus] };
                            if (nP[0] * nC[0] + nP[1] * nC[1] + nP[2] * nC[2] < 0.0f) {
                                nP[0] = -nP[0]; nP[1] = -nP[1]; nP[2] = -nP[2];
                            }
                            if (nM[0] * nC[0] + nM[1] * nC[1] + nM[2] * nC[2] < 0.0f) {
                                nM[0] = -nM[0]; nM[1] = -nM[1]; nM[2] = -nM[2];
                            }
                            out[0] = (nP[0] - nM[0]) * scale;
                            out[1] = (nP[1] - nM[1]) * scale;
                            out[2] = (nP[2] - nM[2]) * scale;
                            };

                        if (i > 0 && i < nx - 1) alignedDiff(idx + 1, idx - 1, inv2dx, &g[0]);
                        else if (i == 0) alignedDiff(idx + 1, idx, invdx, &g[0]);
                        else alignedDiff(idx, idx - 1, invdx, &g[0]);

                        if (j > 0 && j < ny - 1) alignedDiff(idx + strideY, idx - strideY, inv2dy, &g[3]);
                        else if (j == 0) alignedDiff(idx + strideY, idx, invdy, &g[3]);
                        else alignedDiff(idx, idx - strideY, invdy, &g[3]);

                        if (k > 0 && k < nz - 1) alignedDiff(idx + strideZ, idx - strideZ, inv2dz, &g[6]);
                        else if (k == 0) alignedDiff(idx + strideZ, idx, invdz, &g[6]);
                        else alignedDiff(idx, idx - strideZ, invdz, &g[6]);
                    }
                }
            }

            // Evaluate Delta_1 and Delta_2 at node idx using a reference frame (uRef, b1Ref, b2Ref)
            auto evalDeltaAtNode = [&](
                size_t idx,
                const std::array<float, 3>& uRef,
                const std::array<float, 3>& b1Ref,
                float& outD1,
                float& outD2,
                float& outAlign
                ) {
                    float nxVal = nxPtr[idx], nyVal = nyPtr[idx], nzVal = nzPtr[idx];
                    float norm = std::sqrt(nxVal * nxVal + nyVal * nyVal + nzVal * nzVal);
                    if (norm < 1e-8f) {
                        outD1 = 0.0f; outD2 = 0.0f; outAlign = -1.0f;
                        return;
                    }
                    float invN = 1.0f / norm;
                    float n[3] = { nxVal * invN, nyVal * invN, nzVal * invN };

                    float dotRef = n[0] * uRef[0] + n[1] * uRef[1] + n[2] * uRef[2];
                    if (dotRef < 0.0f) {
                        n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
                        dotRef = -dotRef;
                    }
                    outAlign = dotRef;

                    // Parallel-transport b1Ref onto the plane orthogonal to n
                    float proj = b1Ref[0] * n[0] + b1Ref[1] * n[1] + b1Ref[2] * n[2];
                    float d1[3] = { b1Ref[0] - proj * n[0], b1Ref[1] - proj * n[1], b1Ref[2] - proj * n[2] };
                    float d1Len = std::sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
                    if (d1Len < 1e-6f) {
                        std::array<float, 3> nArr = { n[0], n[1], n[2] }, tmp1, tmp2;
                        BuildOrthonormalBasis(nArr, tmp1, tmp2);
                        d1[0] = tmp1[0]; d1[1] = tmp1[1]; d1[2] = tmp1[2];
                    }
                    else {
                        d1[0] /= d1Len; d1[1] /= d1Len; d1[2] /= d1Len;
                    }

                    // d2 = n x d1 completes the right-handed triad {d1, d2, n}
                    float d2[3] = {
                        n[1] * d1[2] - n[2] * d1[1],
                        n[2] * d1[0] - n[0] * d1[2],
                        n[0] * d1[1] - n[1] * d1[0]
                    };

                    const float* g = &gradField[idx * 9];
                    // Directional derivatives (d_alpha . grad) n_b for alpha in {1, 2}
                    float dn1[3], dn2[3];
                    for (int b = 0; b < 3; ++b) {
                        dn1[b] = d1[0] * g[0 + b] + d1[1] * g[3 + b] + d1[2] * g[6 + b];
                        dn2[b] = d2[0] * g[0 + b] + d2[1] * g[3 + b] + d2[2] * g[6 + b];
                    }

                    float d1_dn1 = d1[0] * dn1[0] + d1[1] * dn1[1] + d1[2] * dn1[2];
                    float d2_dn2 = d2[0] * dn2[0] + d2[1] * dn2[1] + d2[2] * dn2[2];
                    float d1_dn2 = d1[0] * dn2[0] + d1[1] * dn2[1] + d1[2] * dn2[2];
                    float d2_dn1 = d2[0] * dn1[0] + d2[1] * dn1[1] + d2[2] * dn1[2];

                    outD1 = 0.5f * (d1_dn1 - d2_dn2);
                    outD2 = 0.5f * (d1_dn2 + d2_dn1);
                };

            auto intersectTri = [&](
                const float pA[3], float f1A, float f2A, float alA,
                const float pB[3], float f1B, float f2B, float alB,
                const float pC[3], float f1C, float f2C, float alC,
                float outPt[3]
                ) -> bool {
                    if ((f1A > 0.0f && f1B > 0.0f && f1C > 0.0f) || (f1A < 0.0f && f1B < 0.0f && f1C < 0.0f) ||
                        (f2A > 0.0f && f2B > 0.0f && f2C > 0.0f) || (f2A < 0.0f && f2B < 0.0f && f2C < 0.0f)) {
                        return false;
                    }
                    if (alA < alignmentThreshold || alB < alignmentThreshold || alC < alignmentThreshold) {
                        return false;
                    }

                    float d1B = f1B - f1A, d1C = f1C - f1A;
                    float d2B = f2B - f2A, d2C = f2C - f2A;
                    float det = d1B * d2C - d1C * d2B;
                    if (std::abs(det) < 1e-12f) return false;

                    float invDet = 1.0f / det;
                    float uB = (-f1A * d2C + f2A * d1C) * invDet;
                    float vC = (-d1B * f2A + d2B * f1A) * invDet;
                    float wA = 1.0f - uB - vC;

                    const float eps = -1e-5f;
                    if (uB < eps || vC < eps || wA < eps) return false;

                    outPt[0] = wA * pA[0] + uB * pB[0] + vC * pC[0];
                    outPt[1] = wA * pA[1] + uB * pB[1] + vC * pC[1];
                    outPt[2] = wA * pA[2] + uB * pB[2] + vC * pC[2];
                    return true;
                };

            static const int faceQuads[6][4] = {
                { 0, 2, 6, 4 },
                { 1, 3, 7, 5 },
                { 0, 1, 5, 4 },
                { 2, 3, 7, 6 },
                { 0, 1, 3, 2 },
                { 4, 5, 7, 6 }
            };

            for (int k = 0; k < nz - 1; ++k) {
                for (int j = 0; j < ny - 1; ++j) {
                    size_t baseRow = strideY * j + strideZ * k;
                    for (int i = 0; i < nx - 1; ++i) {
                        size_t c000 = baseRow + i;
                        size_t cIdx[8] = {
                            c000,
                            c000 + 1,
                            c000 + strideY,
                            c000 + strideY + 1,
                            c000 + strideZ,
                            c000 + strideZ + 1,
                            c000 + strideZ + strideY,
                            c000 + strideZ + strideY + 1
                        };

                        // Quick cell-level rejection in the frame of corner 0
                        std::array<float, 3> n0 = { nxPtr[c000], nyPtr[c000], nzPtr[c000] };
                        float norm0 = std::sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
                        if (norm0 < 1e-8f) continue;
                        std::array<float, 3> uCell = { n0[0] / norm0, n0[1] / norm0, n0[2] / norm0 };
                        std::array<float, 3> b1Cell, b2Cell;
                        BuildOrthonormalBasis(uCell, b1Cell, b2Cell);

                        int maskPos1 = 0, maskNeg1 = 0, maskPos2 = 0, maskNeg2 = 0;
                        bool validCell = true;
                        for (int c = 0; c < 8; ++c) {
                            float d1, d2, al;
                            evalDeltaAtNode(cIdx[c], uCell, b1Cell, d1, d2, al);
                            if (al < alignmentThreshold) { validCell = false; break; }
                            if (d1 >= 0.0f) maskPos1 = 1; else maskNeg1 = 1;
                            if (d2 >= 0.0f) maskPos2 = 1; else maskNeg2 = 1;
                        }
                        if (!validCell || !maskPos1 || !maskNeg1 || !maskPos2 || !maskNeg2) {
                            continue;
                        }

                        float x0 = ox + i * dx, x1 = x0 + dx;
                        float y0 = oy + j * dy, y1 = y0 + dy;
                        float z0 = oz + k * dz, z1 = z0 + dz;

                        const float pos[8][3] = {
                            { x0, y0, z0 },{ x1, y0, z0 },{ x0, y1, z0 },{ x1, y1, z0 },
                            { x0, y0, z1 },{ x1, y0, z1 },{ x0, y1, z1 },{ x1, y1, z1 }
                        };

                        float hitPts[6][3];
                        int hitCount = 0;

                        for (int f = 0; f < 6; ++f) {
                            int v0 = faceQuads[f][0], v1 = faceQuads[f][1], v2 = faceQuads[f][2], v3 = faceQuads[f][3];

                            // Anchor the face frame to v0 so adjacent cells share identical face values
                            size_t refIdx = cIdx[v0];
                            std::array<float, 3> nFace = { nxPtr[refIdx], nyPtr[refIdx], nzPtr[refIdx] };
                            float normF = std::sqrt(nFace[0] * nFace[0] + nFace[1] * nFace[1] + nFace[2] * nFace[2]);
                            if (normF < 1e-8f) continue;
                            std::array<float, 3> uFace = { nFace[0] / normF, nFace[1] / normF, nFace[2] / normF };
                            std::array<float, 3> b1Face, b2Face;
                            BuildOrthonormalBasis(uFace, b1Face, b2Face);

                            int fv[4] = { v0, v1, v2, v3 };
                            float f1[4], f2[4], al[4];
                            for (int m = 0; m < 4; ++m) {
                                evalDeltaAtNode(cIdx[fv[m]], uFace, b1Face, f1[m], f2[m], al[m]);
                            }

                            float pt[3];
                            if (intersectTri(pos[v0], f1[0], f2[0], al[0],
                                pos[v1], f1[1], f2[1], al[1],
                                pos[v2], f1[2], f2[2], al[2], pt) ||
                                intersectTri(pos[v0], f1[0], f2[0], al[0],
                                    pos[v2], f1[2], f2[2], al[2],
                                    pos[v3], f1[3], f2[3], al[3], pt)) {

                                bool duplicate = false;
                                for (int h = 0; h < hitCount; ++h) {
                                    float ddx = pt[0] - hitPts[h][0];
                                    float ddy = pt[1] - hitPts[h][1];
                                    float ddz = pt[2] - hitPts[h][2];
                                    if (ddx * ddx + ddy * ddy + ddz * ddz < 1e-10f) {
                                        duplicate = true;
                                        break;
                                    }
                                }
                                if (!duplicate && hitCount < 6) {
                                    hitPts[hitCount][0] = pt[0];
                                    hitPts[hitCount][1] = pt[1];
                                    hitPts[hitCount][2] = pt[2];
                                    ++hitCount;
                                }
                            }
                        }

                        if (hitCount >= 2) {
                            unsigned int baseIdx = static_cast<unsigned int>(result.vertices.size() / 3);
                            result.vertices.push_back(hitPts[0][0]);
                            result.vertices.push_back(hitPts[0][1]);
                            result.vertices.push_back(hitPts[0][2]);
                            result.vertices.push_back(hitPts[1][0]);
                            result.vertices.push_back(hitPts[1][1]);
                            result.vertices.push_back(hitPts[1][2]);
                            result.lineIndices.push_back(baseIdx);
                            result.lineIndices.push_back(baseIdx + 1);
                        }
                    }
                }
            }

            return result;
        }

        MeshProcessor::ExtractedMesh MeshProcessor::ComputeUmbilicTube(
            const float* directorField,
            float tubeRadius,
            int tubeSides,
            float alignmentThreshold
        ) {
            ExtractedLines lines = ComputeUmbilicLines(directorField, alignmentThreshold);
            return BuildTubeFromExtractedLines(lines, tubeRadius, tubeSides);
        }

        MeshProcessor::ExtractedMesh MeshProcessor::ExtractPolyData(void* polyDataPtr, bool flipOrientation) const {
            auto* mesh = static_cast<vtkPolyData*>(polyDataPtr);
            ExtractedMesh result;

            if (!mesh || mesh->GetNumberOfPoints() == 0 || mesh->GetNumberOfPolys() == 0) {
                return result;
            }

            vtkIdType numPoints = mesh->GetNumberOfPoints();
            result.vertices.resize(numPoints);

            vtkDataArray* ptsData = mesh->GetPoints()->GetData();
            vtkDataArray* normData = mesh->GetPointData()->GetNormals();
            float normSign = flipOrientation ? -1.0f : 1.0f;

            auto* ptsFloat = vtkFloatArray::FastDownCast(ptsData);
            auto* normFloat = vtkFloatArray::FastDownCast(normData);

            if (ptsFloat && normFloat) {
                const float* pPtr = ptsFloat->GetPointer(0);
                const float* nPtr = normFloat->GetPointer(0);
                for (vtkIdType i = 0; i < numPoints; ++i) {
                    auto& v = result.vertices[i];
                    v.position[0] = pPtr[3 * i + 0];
                    v.position[1] = pPtr[3 * i + 1];
                    v.position[2] = pPtr[3 * i + 2];
                    v.normal[0] = normSign * nPtr[3 * i + 0];
                    v.normal[1] = normSign * nPtr[3 * i + 1];
                    v.normal[2] = normSign * nPtr[3 * i + 2];
                }
            }
            else if (ptsFloat) {
                const float* pPtr = ptsFloat->GetPointer(0);
                for (vtkIdType i = 0; i < numPoints; ++i) {
                    auto& v = result.vertices[i];
                    v.position[0] = pPtr[3 * i + 0];
                    v.position[1] = pPtr[3 * i + 1];
                    v.position[2] = pPtr[3 * i + 2];
                    if (normData) {
                        double n[3];
                        normData->GetTuple(i, n);
                        v.normal[0] = normSign * static_cast<float>(n[0]);
                        v.normal[1] = normSign * static_cast<float>(n[1]);
                        v.normal[2] = normSign * static_cast<float>(n[2]);
                    }
                }
            }
            else {
                for (vtkIdType i = 0; i < numPoints; ++i) {
                    double pt[3];
                    mesh->GetPoint(i, pt);
                    auto& v = result.vertices[i];
                    v.position[0] = static_cast<float>(pt[0]);
                    v.position[1] = static_cast<float>(pt[1]);
                    v.position[2] = static_cast<float>(pt[2]);
                    if (normData) {
                        double n[3];
                        normData->GetTuple(i, n);
                        v.normal[0] = normSign * static_cast<float>(n[0]);
                        v.normal[1] = normSign * static_cast<float>(n[1]);
                        v.normal[2] = normSign * static_cast<float>(n[2]);
                    }
                }
            }

            vtkCellArray* polys = mesh->GetPolys();
            result.indices.reserve(polys->GetNumberOfCells() * 3);

            auto cellIter = vtk::TakeSmartPointer(polys->NewIterator());
            for (cellIter->GoToFirstCell(); !cellIter->IsDoneWithTraversal(); cellIter->GoToNextCell()) {
                vtkIdType npts = 0;
                const vtkIdType* pts = nullptr;
                cellIter->GetCurrentCell(npts, pts);
                if (npts >= 3) {
                    unsigned int i0 = static_cast<unsigned int>(pts[0]);
                    unsigned int i1 = static_cast<unsigned int>(pts[1]);
                    unsigned int i2 = static_cast<unsigned int>(pts[2]);
                    if (i0 != i1 && i1 != i2 && i2 != i0) {
                        result.indices.push_back(i0);
                        if (flipOrientation) {
                            result.indices.push_back(i2);
                            result.indices.push_back(i1);
                        }
                        else {
                            result.indices.push_back(i1);
                            result.indices.push_back(i2);
                        }
                    }
                }
            }

            return result;
        }

    }
}