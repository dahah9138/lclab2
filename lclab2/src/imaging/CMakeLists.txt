add_library(lclab2_imaging POM.h POM.cpp RungeSphere.h BMP.h ImageSeries.h ImageSeries.cpp)
target_include_directories(lclab2_imaging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lclab2_imaging PUBLIC lclab2_graphics lclab2_math)
set_property(TARGET lclab2_imaging PROPERTY POSITION_INDEPENDENT_CODE 1)