#ifndef CPUKNN_H
#define CPUKNN_H

#include "Metric.h"
#include <stdio.h>
#include <cstddef>
#include <algorithm>
#include <iostream>
#include <vector>
// Gives std::thread::hardware_concurrency
#include <ppl.h>

/*
 * For each input query point, locates the k-NN (indexes and distances) among the reference points.
 *
 * @param ref        refence points
 * @param ref_nb     number of reference points
 * @param query      query points
 * @param query_nb   number of query points
 * @param dim        dimension of points
 * @param k          number of neighbors to consider
 * @param knn_dist   output array containing the query_nb x k distances
 * @param knn_index  output array containing the query_nb x k indexes
 */

namespace LC { namespace Algorithm {
    using namespace LC::Math;

    namespace CPUKNN {
        // Used for multithreading
        static const std::size_t kNumThreads = 4;// std::thread::hardware_concurrency();
        static int numQueriesProcessed;
        static std::mutex queryLock;
    }

    // Forward declaration of helpers
    template <typename T>
    void modified_insertion_sort(T* dist, int* index, int length, int k);
    template <typename T>
    T compute_distance(const T* ref,
        std::size_t           ref_nb,
        const std::size_t* query,
        const Metric<T>& metric,
        std::size_t           ref_index,
        std::size_t           query_index);
    // ==============================================================

    /**
        * Computes the Euclidean distance between a reference point and a query point.
        *
        * @param ref          refence points
        * @param ref_nb       number of reference points
        * @param query        query points
        * @param query_nb     number of query points
        * @param dim          dimension of points
        * @param ref_index    index to the reference point to consider
        * @param query_index  index to the query point to consider
        * @return computed distance
        */
    template <typename T>
    T compute_distance(const T* ref,
        std::size_t           ref_nb,
        const std::size_t* query,
        const Metric<T>& metric,
        std::size_t           ref_index,
        std::size_t           query_index) {
        std::array<T, 3> r1, r2;

        for (size_t d = 0; d < 3; d++)
        {
            r1[d] = ref[d * ref_nb + ref_index];
            r2[d] = ref[d * ref_nb + query[query_index]];
        }

        return metric.distance(r1, r2);
    }

    /**
        * Gathers at the beginning of the `dist` array the k smallest values and their
        * respective index (in the initial array) in the `index` array. After this call,
        * only the k-smallest distances are available. All other distances might be lost.
        *
        * Since we only need to locate the k smallest distances, sorting the entire array
        * would not be very efficient if k is relatively small. Instead, we perform a
        * simple insertion sort by eventually inserting a given distance in the first
        * k values.
        *
        * @param dist    array containing the `length` distances
        * @param index   array containing the index of the k smallest distances
        * @param length  total number of distances
        * @param k       number of smallest distances to locate
        */
    template <typename T>
    void modified_insertion_sort(T* dist, int* index, int length, int k) {

        // Initialise the first index
        index[0] = 0;

        // Go through all points
        for (int i = 1; i < length; ++i) {

            // Store current distance and associated index
            T curr_dist = dist[i];
            int   curr_index = i;

            // Skip the current value if its index is >= k and if it's higher the k-th slready sorted mallest value
            if (i >= k && curr_dist >= dist[k - 1]) {
                continue;
            }

            // Shift values (and indexes) higher that the current distance to the right
            int j = std::min(i, k - 1);
            while (j > 0 && dist[j - 1] > curr_dist) {
                dist[j] = dist[j - 1];
                index[j] = index[j - 1];
                --j;
            }

            // Write the current distance and index at their position
            dist[j] = curr_dist;
            index[j] = curr_index;
        }
    }

    template <typename T>
    static void knnThread(int start, int end, std::array<bool, 3> B,
        std::size_t k, const T* ref, std::size_t ref_nb, const std::size_t* query, std::size_t query_nb, 
        const Metric<T>& metric, T* knn_dist, std::size_t* knn_index) {

        T* dist = new T[ref_nb];
        int* index = new int[ref_nb];
        size_t refreshRate = 0.001 * query_nb;

        // Lock and check
        CPUKNN::queryLock.lock();
        if (!dist || !index) {
            printf("Memory allocation error\n");
            return;
        }
        CPUKNN::queryLock.unlock();

        // If true then delete after check
        if (!dist || !index) {
            delete[] dist;
            delete[] index;

            return;
        }

        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < ref_nb; ++j) {
                dist[j] = compute_distance<T>(ref, ref_nb, query, metric, j, i);
                index[j] = j;
            }



            // Sort distances / indexes
            // Consider replacing with a gpu parallelizable algorithm
            modified_insertion_sort<T>(dist, index, ref_nb, k);

            // Copy k smallest distances and their associated index
            // No data race, because every thread writes to different pieces of knn_dist, knn_index

            CPUKNN::queryLock.lock();
            for (size_t j = 0; j < k; ++j) {
                if (knn_dist)
                    knn_dist[j * query_nb + i] = dist[j];


                if (knn_index)
                    knn_index[j * query_nb + i] = index[j];

            }

            ++CPUKNN::numQueriesProcessed;

            CPUKNN::queryLock.unlock();

        }

        delete[] dist;
        delete[] index;
    }

    /*
        * For each input query point, locates the k-NN (indexes and distances) among the reference points.
        *
        * @param ref        refence points
        * @param ref_nb     number of reference points
        * @param query      query points
        * @param query_nb   number of query points
        * @param dim        dimension of points
        * @param k          number of neighbors to consider
        * @param knn_dist   output array containing the query_nb x k distances
        * @param knn_index  output array containing the query_nb x k indexes
        */
    template <typename T>
    bool knn_c(const T* ref,
        std::size_t           ref_nb,
        const std::size_t* query,
        std::size_t           query_nb,
        const Metric<T>& metric,
        std::size_t           k,
        T* knn_dist,
        std::size_t* knn_index) {

        std::array<bool, 3> B = metric.Bcs;

        CPUKNN::numQueriesProcessed = 0;

        std::vector<std::thread> threads;
        threads.reserve(CPUKNN::kNumThreads);

        int queriesPerThread = query_nb / CPUKNN::kNumThreads;

        for (int i = 0; i < CPUKNN::kNumThreads; i++) {
            int start = i * queriesPerThread;
            int end = (i == CPUKNN::kNumThreads - 1) ? query_nb : start + queriesPerThread;
            threads.push_back(std::thread(knnThread<T>, start, end, B, k, ref, ref_nb, query, query_nb, std::ref(metric), knn_dist, knn_index));
        }
        for (std::thread& t : threads) t.join();

        return true;

    }

}}

#endif