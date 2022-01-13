#ifndef POM_H
#define POM_H

#include "core.h"
#include "logger.h"
#include "scalar.h"
#include <Eigen/Dense>
#include <complex>

namespace LC { namespace Imaging { namespace UniformGrid {

typedef std::size_t (*t4_2indFunc)(int, int, int, int, const std::array<std::size_t, 3>&);
typedef std::size_t (*t2_2indFunc)(int, int, const std::size_t&);
typedef void (*ColorDataFunc)(void *, const std::array<float, 4>&, std::size_t);

struct POM {
    enum class Waveplate {
		None = 0,
		Full530nm = 1,
        Quarter530nm = 2
	};

    void Compute(scalar *nn, const std::array<int, 3> &voxels, void *CData, ColorDataFunc colorFunc, const float &alpha = 0.5f);

    // Default matlab indexing where slices is of the form
    /*
        { slices[0] = voxels[0],
          slices[1] = voxels[0] * voxels[1],
          slices[2] = voxels[0] * voxels[1] * voxels[2] }
        i.e. a cumulative product of the voxels array
    */

    t4_2indFunc dir2ind = [](int i, int j, int k, int l, const std::array<std::size_t, 3>& slices){
        std::size_t idx = slices[2] * l + slices[1] * k + slices[0] * j + i;
        return idx;
    };

    t2_2indFunc cross2ind = [](int i, int j, const std::size_t& slice) {
        std::size_t idx = j * slice + i;
        return idx;
    };

    Waveplate waveplate = Waveplate::None;
    bool polarizers = true;
    float polarizerAngle = 90.0f;
    double thickness = 0.0;
    double dz = 0.0;
    float gamma = 1.0f;

    float z_scan_ratio = 1.0f;
    int additional_layers = 0;
    int dop = 0;

    double n0 = 0.0;
    double ne = 0.0;

    std::array<float, 3> lightRGB = { 650.0, 550.0, 450.0 };
    std::array<float, 3> intensity = { 1.0, 0.6, 0.2 };

    // Dimension order
    std::array<int, 3> indexOrder = { 0, 1, 2 };



};



}}}




#endif