find_package(ImGui REQUIRED Sources)

set(IMPLOT_HEADERS ${IMPLOT_DIR}/implot.h ${IMPLOT_DIR}/implot_internal.h)
set(IMPLOT_SRC     ${IMPLOT_DIR}/implot.cpp ${IMPLOT_DIR}/implot_items.cpp ${IMPLOT_DIR}/implot_demo.cpp
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


add_library(ImPlot::ImPlot INTERFACE IMPORTED)
set_property(TARGET ImPlot::ImPlot APPEND PROPERTY
	INTERFACE_INCLUDE_DIRECTORIES ${IMPLOT_DIR})
add_library(ImPlot::Sources INTERFACE IMPORTED)
set_property(TARGET ImPlot::Sources APPEND PROPERTY
	INTERFACE_SOURCES ${IMPLOT_SRC})
set_property(TARGET ImPlot::Sources APPEND PROPERTY
	INTERFACE_LINK_LIBRARIES ImPlot::ImPlot)



add_library(LCLAB2ImGuiIntegration SHARED
ImContext.h ImContext.cpp
${ImGuiIntegration_HEADERS} ${ImGuiIntegration_SRCS}
)

target_compile_definitions(LCLAB2ImGuiIntegration PUBLIC IMGUI_USER_CONFIG="${CMAKE_CURRENT_SOURCE_DIR}/visibility.h" PRIVATE LC_BUILD_DLL)
target_link_libraries(LCLAB2ImGuiIntegration
    PUBLIC
        Magnum::GL
        Magnum::Shaders
        ImGui::ImGui
		ImPlot::ImPlot
    PRIVATE
        ImGui::Sources
		ImPlot::Sources)
set_property(TARGET LCLAB2ImGuiIntegration PROPERTY POSITION_INDEPENDENT_CODE 1)
target_include_directories(LCLAB2ImGuiIntegration PUBLIC ${IMGUI_DIR} ${IMPLOT_DIR} ${CMAKE_SOURCE_DIR}/external ${CMAKE_SOURCE_DIR}/external/SDL2-2.0.10/include
															${BASE_DIR}/base)