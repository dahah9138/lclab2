#ifndef FILM_QUADRILATERAL_H
#define FILM_QUADRILATERAL_H

#include "Line.h"

struct Quadrilateral {
    std::array<Vector3, 4> vertices;
    std::array<Line, 4> lines;
    Vector3 lineColor{0.f, 1.f, 0.f}; // Green
    Quadrilateral(const std::array<Vector3, 4>& verts, int *program) : vertices(verts) { Init(*program); }
    Quadrilateral(const std::array<Vector3, 4>& verts, const std::array<int, 4> &ids, int *program) {
        for (int i = 0; i < 4; i++)
            vertices[i] = verts[ids[i]];

        Init(*program);
    }
    void Init(int &program) {
        for (int i = 0; i < 4; i++) {
            int ip1 = (i + 1) % 4;
            lines[i].Create(vertices[i], vertices[ip1], program);
        }
    }

    void Draw(const Matrix4& viewMatrix, const Matrix4& projectionMatrix) {
        Matrix4 mvp = projectionMatrix * viewMatrix;
        for (auto& L : lines)
            L.draw(mvp, lineColor);
    }

    void Draw(Matrix4& mvp) {
        for (auto& L : lines)
            L.draw(mvp, lineColor);
    }
};


#endif