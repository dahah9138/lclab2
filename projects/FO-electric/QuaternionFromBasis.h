#ifndef QUATERNION_FROM_BASIS_H
#define QUATERNION_FROM_BASIS_H

template <typename Ty>
Eigen::Quaternion<Ty> fromBasis(Eigen::Matrix<Ty,3,1> a, Eigen::Matrix<Ty, 3, 1> b, Eigen::Matrix<Ty, 3, 1> c) {
    Ty T = a.x() + b.y() + c.z();
    Ty X, Y, Z, W;
    float s;
    if (T > 0) {
        Ty s = sqrt(T + 1) * 2;
        X = (c.y() - b.z()) / s;
        Y = (a.z() - c.x()) / s;
        Z = (b.x() - a.y()) / s;
        W = 0.25 * s;
    }
    else if (a.x() > b.y() && a.x() > c.z()) {
        s = sqrt(1 + a.x() - b.y() - c.z()) * 2;
        X = 0.25 * s;
        Y = (b.x() + a.y()) / s;
        Z = (a.z() + c.x()) / s;
        W = (c.y() - b.z()) / s;
    }
    else if (b.y() > c.z()) {
        s = sqrt(1 + b.y() - a.x() - c.z()) * 2;
        X = (b.x() + a.y()) / s;
        Y = 0.25 * s;
        Z = (c.y() + b.z()) / s;
        W = (b.z() - c.y()) / s;
    }
    else {
        s = sqrt(1 + c.z() - a.x() - b.y()) * 2;
        X = (a.z() + c.x()) / s;
        Y = (c.y() + b.z()) / s;
        Z = 0.25 * s;
        W = (b.x() - a.y()) / s;
    }
    return Eigen::Quaternion<Ty>{W, X, Y, Z};
}

#endif