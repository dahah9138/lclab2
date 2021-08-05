#ifndef LCLAB2_H
#define LCLAB2_H

// Base
#ifndef LC_CONSOLE_APP
	#include "src/base/Application.h"
#else
	#include "src/base/ConsoleApplication.h"
#endif

#include "src/base/logger.h"
#include "src/base/Header.h"
#include "src/base/Arguments.h"

// LC
#include "src/solver/Solver.h"
#include "src/solver/frankoseen/FOAssets.h"
#include "src/solver/frankoseen/FOFDSolver.h"
#include "src/solver/frankoseen/RBFFDSolver.h"

// Utilities
#include "src/utility/range_pair.h"
#include "src/utility/searchlist.h"

// Math
#include "src/math/vec3.h"
#include "src/math/Choose.h"
#include "src/math/rng.h"
#include "src/math/subset.h"
#include "src/math/rbf.h"
#include "src/math/poly_spline.h"
#include "src/math/powi.h"
#include "src/math/Metric.h"
#include "src/math/AdvancingFront.h"
#include "src/math/StencilWeights.h"

// Algorithms
#include "src/algorithm/cpuknn.h"

#ifndef LC_CONSOLE_APP
	// Graphics
	#include "src/graphics/ArcBall.h"
	#include "src/graphics/SphereArray.h"
	#include "src/graphics/EllipsoidArray.h"
	#include "src/graphics/Sheet.h"
	#include "src/graphics/NormalSheet.h"
	#include "src/graphics/Torus.h"
	#include "src/graphics/NormalTorus.h"
	#include "src/graphics/DynamicColorSheet.h"
	#include "src/graphics/TransparentDrawable.h"

	//implot impl
	#include "src/implementation/ImContext.h"

	// Imaging
	#include "src/imaging/POM.h"
	#include "src/imaging/RungeSphere.h"
#endif

// CUDA
#ifdef LCLAB2_CUDA_AVAIL
    #include "src/cuda/CudaContext.h"
#endif

// ENTRY POINT
#include "src/base/entrypoint.h"

#endif