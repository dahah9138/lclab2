#ifndef LCLAB2_H
#define LCLAB2_H

// Base
#ifndef LC_CONSOLE_APP
	#include "src/base/Application.h"
#else
	#include "src/base/ConsoleApplication.h"
#endif

#include "src/base/scalar.h"

#include "src/base/logger.h"
#include "src/base/Header.h"
#include "src/base/Arguments.h"

// LC
#include "src/solver/Solver.h"
#include "src/solver/frankoseen/FOAssets.h"
#include "src/solver/frankoseen/FOFDSolver.h"
#include "src/solver/frankoseen/RBFFDSolver.h"
#include "src/solver/film/FilmSolver.h"

#include "src/solver/qtensor/QTensorAssets.h"
#include "src/solver/qtensor/QTensorSolver.h"

// Utilities
#include "src/utility/range_pair.h"
#include "src/utility/searchlist.h"

// Math
#include "src/math/Choose.h"
#include "src/math/rng.h"
#include "src/math/subset.h"
#include "src/math/rbf.h"
#include "src/math/poly_spline.h"
#include "src/math/powi.h"
#include "src/math/Metric.h"
#include "src/math/LinearInterpolate.h"
#include "src/math/LinearInterpolator.h"
#include "src/math/Isosurface.h"
#include "src/math/ExtendedMC/MarchingCubes.h"
#include "src/math/AdvancingFront.h"
#include "src/math/StencilWeights.h"
#include "src/math/CumulativeTrapIntegral.h"
#include "src/math/HopfCharge.h"
#include "src/math/BaryonDensity.h"
#include "src/math/MaxEigen.h"
#include "src/math/ChiralityTensor.h"
#include "src/math/Derivative.h"
#include "src/math/ChiralityField.h"
#include "src/math/ScalarOrderParameter.h"
#include "src/math/Graph.h"
#include "src/math/Normalize.h"

// Algorithms
#include "src/algorithm/cpuknn.h"
#include "src/algorithm/loop_subdivision/loop_subdivision_impl.h"

// Smoothing algorithm
#include "src/smoothing/smooth_alg.hpp"

#ifndef LC_CONSOLE_APP
	// Graphics
	#include "src/graphics/ArcBall.h"
	#include "src/graphics/SphereArray.h"
	#include "src/graphics/EllipsoidArray.h"
	#include "src/graphics/NematicArray.h"
	#include "src/graphics/Sheet.h"
	#include "src/graphics/NormalSheet.h"
	#include "src/graphics/Torus.h"
	#include "src/graphics/NormalTorus.h"
	#include "src/graphics/TubularSurface.h"
	#include "src/graphics/DynamicColorSheet.h"
	#include "src/graphics/Surface.h"
	#include "src/graphics/TransparentDrawable.h"
	#include "src/graphics/TransparentNormalDrawable.h"
	#include <Magnum/Primitives/Line.h>

	//implot impl
	#include "src/implementation/ImContext.h"
#endif

// Imaging
#include "src/imaging/POM.h"
#include "src/imaging/RungeSphere.h"

#include "src/imaging/BMP.h"
#include "src/imaging/ImageSeries.h"

// CUDA
#ifdef LCLAB2_CUDA_AVAIL
    #include "src/cuda/CudaContext.h"
#endif

// ENTRY POINT
#include "src/base/entrypoint.h"

#endif