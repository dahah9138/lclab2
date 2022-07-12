#ifndef TUBULARSURFACE_H
#define TUBULARSURFACE_H

#include "NormalSheet.h"
#include <eigen/Eigen/Dense>

// Based on
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.104.7190&rep=rep1&type=pdf

namespace LC {

struct TubularSurface: public NormalSheet {

    void Init(const bool & closedTube);


    // Points used to create tubular surface
    std::vector< Eigen::Vector3d > curve;
    std::vector<Magnum::UnsignedInt> indices;
    float radius = 0.25f;
    int relaxIterations = 40;
    float alpha = 1.f;
    bool compile = false;
};

}


#endif
