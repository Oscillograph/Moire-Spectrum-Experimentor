#ifndef BSE_COMMON_STL_H
#define BSE_COMMON_STL_H

// interfaces and algorithms
#include <cstdlib> // exit and other functions
#include <stdio.h> // C I/O
#include <iostream> // C++ I/O
#include <fstream> // file input/output
#include <memory> // smart pointers
#include <utility> // pairs
#include <algorithm> // search functions and others
#include <functional> // not decided yet, but Events system uses it, so why not
#include <initializer_list> // to convert initializer lists to vectors to avoid problems with pointers to vectors
#include <chrono> // time management
#include <math.h> // simple math for graphics, sound and other tasks

// types
#include <string> // text operations, obviously
#include <sstream> // for easier work with strings
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>

#ifdef SAVANNAH_PLATFORM_WINDOWS
	#include <windows.h>
#endif

// vendor libs
/*
#include <glm/glm.hpp> // complex math for graphics
#include <glm/gtc/matrix_transform.hpp> // matrix transform and translate libraries
#include <glm/gtc/type_ptr.hpp> // for value_ptr function
#include <glm/gtx/quaternion.hpp> // for quaternions
#include <glm/gtx/matrix_decompose.hpp> // for glm::decompose

#include <imgui/imgui.h> // in systems/gui/ImGuiLayer.h
#include <imgui/imgui_internal.h> // to use PushMultiItemsWidths function, etc.
#include <imgui/misc/cpp/imgui_stdlib.h> // std::string wrappers to make life easier

#include <yaml_cpp/include/yaml-cpp/yaml.h> // serialization/deserialization library

#include <external/utf8.h> // to make it easier to work with utf8, utf32
*/

#endif
