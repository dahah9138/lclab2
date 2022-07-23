#ifndef LEHMAN_CLUSTER_H
#define LEHMAN_CLUSTER_H

#include <lclab2.h>

/*
    Contains a lehman cross section that can be embedded in a uniform helical bg
*/

typedef unsigned int uint;

struct lehman_cross_section {
    struct Node {
        Node() = default;
        Node(const Eigen::Vector3d& p, const Eigen::Vector3d& d) : position(p), director(d) {}
        Eigen::Vector3d position;
        Eigen::Vector3d director;
    };

    // Fill the Lehman cluster cross section according to resolution provided by cell/(vox-1)
    void fill(std::array<int, 3> vox, const std::array<LC::scalar, 3>& cell, int upsample) {

        // Clear data
        nodes.clear();
        cell_dims.x() = cell[0];
        cell_dims.y() = cell[1];
        cell_dims.z() = cell[2];
        voxels = vox;

        for (int i = 0; i < 3; i++)
            vox[i] *= upsample;

        auto helical = LC::Math::Planar(2, 1);

        LC::scalar dx = cell[0] / (vox[0] - 1);
        LC::scalar dz = cell[2] / (vox[2] - 1);

        // Orientation is initally yhat
        orientation = { 0., 1., 0. };

        // Input:
        // r0: Position of the defect
        // y0: Y plane where the defect is placed
        // sgn: Which side the defect is placed (+ right, - left)
        auto LehmanCluster = [&](Eigen::Vector3d r0, int sgn) {

            // Only works for z = 0 defect location right now
            std::array<LC::scalar, 3> n0 = helical(r0.x(), r0.y(), 0.);

            for (int x = 0; x < vox[0]; x++) {
                for (int z = 0; z < vox[2]; z++) {
                    // Get current position
                    Eigen::Vector3d pos(
                        -cell[0] * 0.5 + x * dx,
                        0.,
                        -cell[2] * 0.5 + z * dz
                    );

                    // Position relative to the point defect
                    Eigen::Vector3d rprime = pos - r0;
                    LC::scalar rprime_len = rprime.norm();

                    // If within half a pitch from the point defect and to the right of the defect
                    if (rprime_len > 0. && rprime_len <= 0.5 && sgn * rprime.x() > 0.) {
                        Eigen::Quaterniond yhat(0., n0[0], n0[1], n0[2]); // Initial director orientation at center of defect
                        // Angle of position relative to point defect in the defect plane
                        LC::scalar phi = atan2(rprime.z(), rprime.x());

                        // Radial rotation quaternion
                        LC::scalar rprime_len = rprime.norm();
                        LC::scalar theta = 2. * M_PI * rprime_len;
                        LC::scalar ct, st;
                        ct = cos(0.5 * theta);
                        st = sin(0.5 * theta);
                        Eigen::Quaterniond rot_quat;
                        rot_quat.w() = ct;
                        rot_quat.x() = st * rprime.x() / rprime_len;
                        rot_quat.y() = st * rprime.y() / rprime_len;
                        rot_quat.z() = st * rprime.z() / rprime_len;

                        // Apply the quaternion to yhat
                        Eigen::Quaterniond n = rot_quat * yhat * rot_quat.conjugate();
                        Eigen::Vector3d nn;
                        nn(0) = n.x();
                        nn(1) = n.y();
                        nn(2) = n.z();
                        nodes.emplace_back(pos, nn);
                    }
                }
            }
        };

        // Create the Lehman clusters starting in the y = 0 plane
        // -------------------------------------------------------
        // Lehman cluster at x = -0.5 curved to the right
        LehmanCluster({ -0.5, 0., 0. }, 1);
        // Lehman cluster at x = 0.5 curved to the left
        LehmanCluster({ 0.5, 0., 0. }, -1);

        reset();
    }

    // Data transformations

    // Translate the Lehman cluster forward by translation_units
    void forward(LC::scalar translation_units) {

        // Need rotation to track orientation change
        LC::scalar ct = cos(total_rotation * 0.5);
        LC::scalar st = sin(total_rotation * 0.5);
        Eigen::Quaterniond quat(ct, 0., 0., st); // cos(theta/2) + sin(theta/2) * zhat
        Eigen::Matrix3d rot = quat.toRotationMatrix();

        // Apply rotation to the plane orientation
        Eigen::Vector3d rot_orientation = rot * orientation;
        center = center + translation_units * rot_orientation;

        for (auto & node : transformed_nodes) {
            node.position = node.position + translation_units * rot_orientation;
        }
    }

    // Translate the Lehman cluster in the z-direction by dist dz
    // path length is the distance traveled in the new orientation of the plane determined by dz
    void rotation(LC::scalar path_len, LC::scalar dz) {
        // Create rotation matrix
        LC::scalar theta = 2 * M_PI * dz;
        total_rotation += theta;
        LC::scalar ct = cos(total_rotation * 0.5);
        LC::scalar st = sin(total_rotation * 0.5);
        Eigen::Quaterniond quat(ct, 0., 0., st); // cos(theta/2) + sin(theta/2) * zhat
        Eigen::Matrix3d rot = quat.toRotationMatrix();

        // Apply rotation to the plane orientation
        Eigen::Vector3d rot_orientation = rot * orientation;

        // Apply quaternion to rotate positions
        for (uint i = 0; i < transformed_nodes.size(); i++) {
            // Can apply rotation immediately to directors
            transformed_nodes[i].director = rot * nodes[i].director;
            // Translate back to center and apply rotation, then translate back 
            // and translate along plane orientation by amount path_len
            transformed_nodes[i].position = rot * nodes[i].position + center + path_len * rot_orientation;
            transformed_nodes[i].position.z() += dz;
        }

        center = center + path_len * rot_orientation + Eigen::Vector3d(0., 0., dz);
    }

    // Copy initial data to transformed data again
    void reset() {
        total_rotation = 0.;
        center = Eigen::Vector3d(0., 0., 0.);
        transformed_nodes.clear();
        transformed_nodes.reserve(nodes.size());
        for (const auto& node : nodes)
            transformed_nodes.emplace_back(node);
    }

    // Extract the index for the transformed position at index id in the data
    uint index(uint id) {
        int i = (transformed_nodes[id].position.x() / cell_dims.x() + 0.5) * (voxels[0] - 1);
        int j = (transformed_nodes[id].position.y() / cell_dims.y() + 0.5) * (voxels[1] - 1);
        int k = (transformed_nodes[id].position.z() / cell_dims.z() + 0.5) * (voxels[2] - 1);
        return i + voxels[0] * j + voxels[0] * voxels[1] * k;
    }
    
    uint size() {
        return nodes.size();
    };

    // Cross section data
    std::vector<Node> nodes;

    Eigen::Vector3d cell_dims, center;
    LC::scalar total_rotation;

    std::array<int, 3> voxels;

    // Orientation of the plane embedding the Lehman cluster cross section
    Eigen::Vector3d orientation;

    // Transformed nodes
    std::vector<Node> transformed_nodes;
};


#endif