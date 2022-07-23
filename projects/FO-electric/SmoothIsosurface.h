#ifndef SMOOTH_ISOSURFACE_H
#define SMOOTH_ISOSURFACE_H

#include <lclab2.h>

// Smooth isosurfaces
void SmoothIsosurface(LC::Math::IsoVertex* verts, unsigned int* indices, unsigned int nVert, unsigned int nInd, int iterations, float smoothingValue, int smoothingType) {
    unsigned int nTriangles = nInd / 3;

    if (iterations < 1)
        return;

    // Extract vertex data
    Smoothing::Mesh mesh;
    mesh._vertices.resize(nVert);
    mesh._normals.resize(nVert);
    mesh._triangles.resize(nTriangles);

    for (int i = 0; i < nVert; i++) {
        mesh._vertices[i].x = verts[i].position[0];
        mesh._vertices[i].y = verts[i].position[1];
        mesh._vertices[i].z = verts[i].position[2];
    }

    for (int i = 0; i < nTriangles; i++) {
        mesh._triangles[i][0] = indices[3 * i];
        mesh._triangles[i][1] = indices[3 * i + 1];
        mesh._triangles[i][2] = indices[3 * i + 2];
    }

    Smoothing::Vertex_to_1st_ring_vertices first_ring;
    Smoothing::Vertex_to_face v_to_face;
    v_to_face.compute(mesh);
    first_ring.compute(mesh, v_to_face);
    
    if (smoothingType == 0)
        mesh._vertices = Smoothing::explicit_laplacian_smoothing(mesh._vertices, first_ring._rings_per_vertex, iterations, smoothingValue);
    else if (smoothingType == 1)
        mesh._vertices = Smoothing::smooth_iterative(mesh._vertices, first_ring._rings_per_vertex, iterations, smoothingValue);
    else if (smoothingType == 2)
        mesh._vertices = Smoothing::implicit_laplacian_smoothing(mesh._vertices, first_ring._rings_per_vertex, iterations, smoothingValue);

    // Copy data back
    for (int i = 0; i < nVert; i++) {
        verts[i].position[0] = mesh._vertices[i].x;
        verts[i].position[1] = mesh._vertices[i].y;
        verts[i].position[2] = mesh._vertices[i].z;
    }

    // Set normals to zero
    for (unsigned int i = 0; i < nVert; i++) {
        verts[i].normal[0] = 0;
        verts[i].normal[1] = 0;
        verts[i].normal[2] = 0;
    }

    // Recompute normals
    for (int i = 0; i < nTriangles; i++) {

        LC::Math::VECTOR3D vec1, vec2, normal;
        unsigned int id0, id1, id2;
        id0 = indices[i * 3];
        id1 = indices[i * 3 + 1];
        id2 = indices[i * 3 + 2];
        vec1[0] = verts[id1].position[0] - verts[id0].position[0];
        vec1[1] = verts[id1].position[1] - verts[id0].position[1];
        vec1[2] = verts[id1].position[2] - verts[id0].position[2];
        vec2[0] = verts[id2].position[0] - verts[id0].position[0];
        vec2[1] = verts[id2].position[1] - verts[id0].position[1];
        vec2[2] = verts[id2].position[2] - verts[id0].position[2];

        LC::Math::CrossProduct(vec1, vec2, normal);

        verts[id0].normal[0] += normal[0];
        verts[id0].normal[1] += normal[1];
        verts[id0].normal[2] += normal[2];
        verts[id1].normal[0] += normal[0];
        verts[id1].normal[1] += normal[1];
        verts[id1].normal[2] += normal[2];
        verts[id2].normal[0] += normal[0];
        verts[id2].normal[1] += normal[1];
        verts[id2].normal[2] += normal[2];
    }

    // Normalize normals
    for (unsigned int i = 0; i < nVert; i++) {
        float len = 0.0f;
        for (int d = 0; d < 3; d++)
            len += verts[i].normal[d] * verts[i].normal[d];

        len = sqrt(len);

        verts[i].normal[0] /= len;
        verts[i].normal[1] /= len;
        verts[i].normal[2] /= len;
    }
}


#endif