#pragma once

// Base
#include "src/base/Application.h"
#include "src/base/logger.h"
#include "src/base/Header.h"

// LC
#include "src/solver/Solver.h"
#include "src/solver/frankoseen/FOAssets.h"
#include "src/solver/frankoseen/FOFDSolver.h"
#include "src/solver/frankoseen/RBFFDSolver.h"
// Math
#include "src/math/vec3.h"
// Graphics
#include "src/graphics/ArcBall.h"
#include "src/graphics/SphereArray.h"
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

// CUDA
#ifdef LCLAB2_CUDA
    #include "src/cuda/CudaContext.h"
#endif

// ENTRY POINT
#include "src/base/entrypoint.h"
