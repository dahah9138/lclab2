#ifndef MATH_FIELDS_H
#define MATH_FIELDS_H

#include "Configuration.h"
#include <Eigen/Dense>

namespace LC { namespace Math {

	using namespace Configuration;

	// Director configurations
	VectorField Uniform(std::array<scalar, 3> n0 = { 0.0, 0.0, 1.0 });
	VectorField Planar(int layers, scalar cellZ);
	VectorField Heliknoton(int Q, std::array<scalar, 3> cell, scalar lambda = 1.0, scalar lim = 1.135, const Eigen::Matrix<scalar, 3, 1>& translation = Eigen::Matrix<scalar, 3, 1>{ 0.0, 0.0, 0.0 }, bool background = true);
	VectorField Hopfion(int Q, std::array<scalar, 3> cell, scalar lambda = 1.0, scalar lim = 2.0, bool background = true);


	// Exclusion radii
	ScalarField UniformRadius(scalar r = 0.1);
	ScalarField LinearSphere(scalar r1, scalar r2, scalar rend, scalar rstart = 0.0);
	ScalarField YLine(scalar r1, scalar r2, scalar xstart, scalar xend);
	ScalarField ZLine(scalar r1, scalar r2, scalar zstart, scalar zend);

	// Voltage functions
	ScalarField VoltageZ(scalar vi, scalar vf, scalar cellZ);
	
	// Active regions
	IsActive ActiveSphere(scalar r);
	IsActive ActiveParallelopiped(std::array<scalar, 3> cell, std::array<bool, 3> bd);
	// Check if in torus created by (r1 - sqrt(x1^2+x2^2))^2 + x3^2 = r^2
	// where x_i = (position - (x,y,z))_i
	IsActive Torus(std::array<scalar, 3> position, scalar r1, scalar r2);
	IsActive RubberTorus(std::array<scalar, 3> position, scalar r1, scalar r2, scalar deformation = 0.0);
	
	
}}


#endif