#include "cpuknn.h"

namespace LC { namespace Algorithm {
    using namespace LC::Math;

    // Templated functions

    template <>
    float compute_distance(const float* ref,
        std::size_t           ref_nb,
        const std::size_t* query,
        std::size_t query_nb,
        const Metric<float>& metric,
        std::size_t           ref_index,
        std::size_t           query_index) {

        std::array<float, 3> r1, r2;

        for (size_t d = 0; d < 3; d++)
        {
            r1[d] = ref[d * ref_nb + ref_index];
            r2[d] = ref[d * ref_nb + query[query_index]];
        }

        return metric.distance(r1, r2);
    }

    template <>
    float compute_distance(const float* ref,
        std::size_t           ref_nb,
        const float* query,
        std::size_t query_nb,
        const Metric<float>& metric,
        std::size_t           ref_index,
        std::size_t           query_index) {

        std::array<float, 3> r1, r2;

        for (size_t d = 0; d < 3; d++)
        {
            r1[d] = ref[d * ref_nb + ref_index];
            r2[d] = query[d * query_nb + query_index];
        }

        return metric.distance(r1, r2);
    }


    template <>
    double compute_distance(const double* ref,
        std::size_t           ref_nb,
        const std::size_t* query,
        std::size_t query_nb,
        const Metric<double>& metric,
        std::size_t           ref_index,
        std::size_t           query_index) {

        std::array<double, 3> r1, r2;

        for (size_t d = 0; d < 3; d++)
        {
            r1[d] = ref[d * ref_nb + ref_index];
            r2[d] = ref[d * ref_nb + query[query_index]];
        }

        return metric.distance(r1, r2);
    }

    template <>
    double compute_distance(const double* ref,
        std::size_t           ref_nb,
        const double* query,
        std::size_t query_nb,
        const Metric<double>& metric,
        std::size_t           ref_index,
        std::size_t           query_index) {

        std::array<double, 3> r1, r2;

        for (size_t d = 0; d < 3; d++)
        {
            r1[d] = ref[d * ref_nb + ref_index];
            r2[d] = query[d * query_nb + query_index];
        }

        return metric.distance(r1, r2);
    }



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

    template<> bool knn_c(const float* ref,
        std::size_t           ref_nb,
        const float* query,
        std::size_t           query_nb,
        const Metric<float>& metric,
        std::size_t           k,
        float* knn_dist,
        std::size_t* knn_index);

    template<> bool knn_c(const double* ref,
        std::size_t           ref_nb,
        const double* query,
        std::size_t           query_nb,
        const Metric<double>& metric,
        std::size_t           k,
        double* knn_dist,
        std::size_t* knn_index);

}}