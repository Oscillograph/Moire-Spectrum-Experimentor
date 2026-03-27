cmake_minimum_required(VERSION 3.28)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(SAVANNAH_IMGUI_DEV "Dear ImGui development library" ON)

# set the project name
project(moiree-spectrum-experimentor
    VERSION 1.0.0
    DESCRIPTION "Experiments with Moiree patterns on LFM-CW spectrums"
    LANGUAGES CXX)

# global strings
set(PROJECT_BASE_SOURCE_DIR ${CMAKE_SOURCE_DIR})
message("set PROJECT_BASE_SOURCE_DIR to ${PROJECT_BASE_SOURCE_DIR}")

set(SAVANNAH_INCLUDES_DIR 
    ${PROJECT_BASE_SOURCE_DIR}/include
    ${PROJECT_BASE_SOURCE_DIR}/external
    ${PROJECT_BASE_SOURCE_DIR}/external/imgui
    ${PROJECT_BASE_SOURCE_DIR}/external/implot
    ${PROJECT_BASE_SOURCE_DIR}/external/glad/include
    ${PROJECT_BASE_SOURCE_DIR}/external/imgui_filedialog/src
    ${PROJECT_BASE_SOURCE_DIR}/
    )
    
# global includes
include_directories(${SAVANNAH_INCLUDES_DIR})

# link external libraries
#add_library(savannah STATIC IMPORTED)
#set_target_properties(savannah PROPERTIES IMPORTED_LOCATION ${PROJECT_BASE_SOURCE_DIR}/bin/linux/libsavannah.a)
#target_include_directories(savannah INTERFACE ${PROJECT_BASE_SOURCE_DIR}/include)

add_library(glfw3 STATIC IMPORTED)
set_target_properties(glfw3 PROPERTIES IMPORTED_LOCATION ${PROJECT_BASE_SOURCE_DIR}/bin/linux/libglfw3.a)
target_include_directories(glfw3 INTERFACE ${PROJECT_BASE_SOURCE_DIR}/external/glfw)

add_library(glad STATIC IMPORTED)
set_target_properties(glad PROPERTIES IMPORTED_LOCATION ${PROJECT_BASE_SOURCE_DIR}/bin/linux/libglad.a)
target_include_directories(glad INTERFACE ${PROJECT_BASE_SOURCE_DIR}/external/glad/include)

#add_library(dear-imgui STATIC IMPORTED)
#if(SAVANNAH_IMGUI_DEV)
#	set_target_properties(dear-imgui PROPERTIES IMPORTED_LOCATION ${PROJECT_BASE_SOURCE_DIR}/bin/linux/libimgui-opengl-dev.a)
#else()
#	set_target_properties(dear-imgui PROPERTIES IMPORTED_LOCATION ${PROJECT_BASE_SOURCE_DIR}/bin/linux/libimgui-opengl.a)
#endif()
#target_include_directories(dear-imgui INTERFACE ${PROJECT_BASE_SOURCE_DIR}/external/imgui)

add_library(yaml-cpp STATIC IMPORTED)
set_target_properties(yaml-cpp PROPERTIES IMPORTED_LOCATION ${PROJECT_BASE_SOURCE_DIR}/bin/linux/libyaml-cpp.a)
target_include_directories(yaml-cpp INTERFACE ${PROJECT_BASE_SOURCE_DIR}/external/yaml_cpp/include)

add_library(fftw STATIC IMPORTED)
set_target_properties(fftw PROPERTIES IMPORTED_LOCATION ${PROJECT_BASE_SOURCE_DIR}/bin/linux/libfftw3.a)
target_include_directories(fftw INTERFACE ${PROJECT_BASE_SOURCE_DIR}/external/fftw)

include(${PROJECT_BASE_SOURCE_DIR}/cmake/build-script-linux.cmake)


# ======================================
# Project-specific
# --------------------------------------
savannah_add_target(
	${PROJECT_NAME}
	EXE
	"${PROJECT_BASE_SOURCE_DIR}/include
	 ${PROJECT_BASE_SOURCE_DIR}"
	"${PROJECT_BASE_SOURCE_DIR}/source 
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/imgui.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/implot/implot.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/implot/implot_items.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/implot/implot_demo.cpp
 	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/misc/cpp/imgui_stdlib.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/imgui_draw.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/imgui_demo.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/imgui_tables.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/imgui_widgets.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/backends/imgui_impl_opengl3.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui/backends/imgui_impl_glfw.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/imgui_filedialog/src/imgui_filedialog.cpp
	 ${PROJECT_BASE_SOURCE_DIR}/external/stb_image.cpp
	 "
	${PROJECT_BASE_SOURCE_DIR}/out
	".."
	"-Wall -Wextra -DSAVANNAH_LOGGER_CLASS"
	)
