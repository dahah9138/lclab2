# Need to compile glad loader
add_library(GLADSTATICLIB ${CMAKE_SOURCE_DIR}/external/glad/glad.c)
# For glad.h
target_include_directories(GLADSTATICLIB PRIVATE ${CMAKE_SOURCE_DIR}/external)
set_target_properties (GLADSTATICLIB PROPERTIES LINKER_LANGUAGE C)

# Compile ImPlot
#add_library(ImPlot ${IMPLOT_DIR}/implot.h ${IMPLOT_DIR}/implot.cpp
#			${IMPLOT_DIR}/implot_internal.h ${IMPLOT_DIR}/implot_items.cpp)
# One of them should have glad...
#target_link_libraries(ImPlot PRIVATE MagnumIntegration::ImGui
#									Magnum::GL
#									Magnum::Magnum
#									Magnum::Shaders)
#target_compile_definitions(ImPlot PRIVATE
#IMGUI_IMPL_OPENGL_LOADER_GLAD
#IMPLOT_API
#)


#if (WIN32)
#	target_include_directories(ImPlot PRIVATE ${CMAKE_SOURCE_DIR}/external/SDL2-2.0.10/include ${IMGUI_DIR} ${CMAKE_SOURCE_DIR}/external)
#else()
#	target_include_directories(ImPlot PRIVATE ${IMGUI_DIR} ${CMAKE_SOURCE_DIR}/external)
#endif()


set(IMGUI_HEADERS
    ${IMGUI_DIR}/imconfig.h
    ${IMGUI_DIR}/imgui_internal.h
    ${IMGUI_DIR}/imgui.h
    ${IMGUI_DIR}/imstb_rectpack.h
    ${IMGUI_DIR}/imstb_textedit.h
    ${IMGUI_DIR}/imstb_truetype.h
    ${IMGUI_DIR}/misc/cpp/imgui_stdlib.h
)

set(IMGUI_SRC
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/misc/cpp/imgui_stdlib.cpp
    ${IMGUI_DIR}/backends/imgui_impl_sdl.cpp
    ${IMGUI_DIR}/backends/imgui_impl_sdl.h
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.h
)

# Copy pasted directly from magnum-plugins...
set(ImGuiIntegration_SRCS
    Context.cpp)

set(ImGuiIntegration_HEADERS
    Context.h
    Context.hpp
    Integration.h
    Widgets.h

    visibility.h)
	
set(IMPLOT_HEADERS ${IMPLOT_DIR}/implot.h ${IMPLOT_DIR}/implot_internal.h)
set(IMPLOT_SRC     ${IMPLOT_DIR}/implot.cpp ${IMPLOT_DIR}/implot_items.cpp #${IMPLOT_DIR}/implot_demo.cpp
)

add_library(LCLAB2ImplementationImLib ImContext.h ImContext.cpp ${ImGuiIntegration_HEADERS} ${ImGuiIntegration_SRCS})
target_sources(LCLAB2ImplementationImLib PRIVATE ${IMPLOT_SRC} ${IMPLOT_HEADERS})
target_sources(LCLAB2ImplementationImLib PRIVATE ${IMGUI_SRC} ${IMGUI_HEADERS})
#target_link_libraries(LCLAB2ImplementationImPlot PUBLIC MagnumIntegration::ImGui)
target_link_libraries(LCLAB2ImplementationImLib PUBLIC
													Magnum::GL
													Magnum::Shaders
												PRIVATE
													GLADSTATICLIB)
target_compile_definitions(LCLAB2ImplementationImLib PRIVATE MAG_BUILD_DLL)

target_include_directories(LCLAB2ImplementationImLib PUBLIC ${IMGUI_DIR} ${IMPLOT_DIR} ${CMAKE_SOURCE_DIR}/external ${CMAKE_SOURCE_DIR}/external/SDL2-2.0.10/include
															${BASE_DIR}/base)