add_library(lclab2_math Choose.h Choose.cpp 
						rng.h rng.cpp 
						subset.h subset.cpp 
						rbf.h 
						poly_spline.h poly_spline.cpp 
						powi.h powi.cpp
						Metric.h Metric.cpp
						AdvancingFront.h AdvancingFront.cpp
						StencilWeights.h StencilWeights.cpp
						LinearInterpolate.h LinearInterpolate.cpp
						LinearInterpolator.h LinearInterpolator.cpp
						Vector.h Vector.cpp
						HopfCharge.h BaryonDensity.h
						CumulativeTrapIntegral.h CumulativeTrapIntegral.cpp
						MaxEigen.h
						ChiralityTensor.h
						Derivative.h
						ChiralityField.h ChiralityField.cpp
						Isosurface.h Isosurface.cpp
						ExtendedMC/defs.h ExtendedMC/LookUpTable.h
						ExtendedMC/MarchingCubes.h ExtendedMC/MarchingCubes.cpp
						ExtendedMC/ply.h ExtendedMC/ply.cpp
						ScalarOrderParameter.h
						Graph.h
						Normalize.h
						Fields.h Fields.cpp)
target_include_directories(lclab2_math PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_SOURCE_DIR}/lclab2/src/base
${CMAKE_SOURCE_DIR}/external/eigen
${CMAKE_SOURCE_DIR}/external/eigen/unsupported)
target_link_libraries(lclab2_math PUBLIC lclab2_utility)
set_property(TARGET lclab2_math PROPERTY POSITION_INDEPENDENT_CODE 1)