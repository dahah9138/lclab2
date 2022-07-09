#ifndef SMOOTHING_SOLVERS_HPP
#define SMOOTHING_SOLVERS_HPP

#include <vector>
#include "mesh.hpp"
#include "vec3.hpp"
#include "vertex_to_1st_ring_vertices.hpp"

// https://github.com/brainexcerpts/laplacian_smoothing_triangle_mesh

namespace Smoothing {
/// @param alpha : strength of the smoothing, with the explicit scheme it should
/// be alpha < 1 and often I recommend alpha = 0.5 when using cotangent weights)... :/
std::vector< Vec3 >
explicit_laplacian_smoothing(const std::vector<Vec3>& in_vertices,
                             const std::vector< std::vector<int> >& edges,
                             int nb_iter,
                             float alpha = 0.5f);

/// Equivalent of the explicit_laplacian_smoothing() but we efficiently
/// iterate over the mesh instead of the matrix representation.
std::vector< Vec3 >
smooth_iterative(const std::vector<Vec3>& in_vertices,
                 const std::vector< std::vector<int> >& edges,
                 int nb_iter,
                 float alpha = 0.5f);

/// @param alpha : strength of the smoothing, with the implicit scheme it can be
/// arbitrary large! :D
std::vector< Vec3 >
implicit_laplacian_smoothing(const std::vector<Vec3>& in_vertices,
                             const std::vector< std::vector<int> >& edges,
                             int nb_iter,
                             float alpha = 100.0f);
}

#endif // SMOOTHING_SOLVERS_HPP
