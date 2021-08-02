#include "cpuknn.h"

namespace LC { namespace Algorithm {
    using namespace LC::Math;

    // Templated functions
    template<> bool knn_c(const float* ref,
        unsigned int           ref_nb,
        const unsigned int* query,
        unsigned int           query_nb,
        const Metric<float>& metric,
        unsigned int           k,
        float* knn_dist,
        unsigned int* knn_index);

    template<> bool knn_c(const double* ref,
        unsigned int           ref_nb,
        const unsigned int* query,
        unsigned int           query_nb,
        const Metric<double>& metric,
        unsigned int           k,
        double* knn_dist,
        unsigned int* knn_index);

}}