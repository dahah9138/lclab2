#include "cpuknn.h"

namespace LC { namespace Algorithm {
    using namespace LC::Math;

    // Templated functions
    template<> bool knn_c(const float* ref,
        std::size_t           ref_nb,
        const std::size_t* query,
        std::size_t           query_nb,
        const Metric<float>& metric,
        std::size_t           k,
        float* knn_dist,
        std::size_t* knn_index);

    template<> bool knn_c(const double* ref,
        std::size_t           ref_nb,
        const std::size_t* query,
        std::size_t           query_nb,
        const Metric<double>& metric,
        std::size_t           k,
        double* knn_dist,
        std::size_t* knn_index);

}}