#ifndef SAVANNAH_H
#define SAVANNAH_H

#include <savannah/core.h>

#include <savannah/main-app.h>

// #include <savannah/core/application/application.h>
// #include <savannah/core/keycodes.h>

// Savannah systems
// - events system
// #include <savannah/systems/events/event.h>
// #include <savannah/systems/events/appEvent.h>
// #include <savannah/systems/events/keyEvent.h>
// #include <savannah/systems/events/mouseEvent.h>

// - scenes
// #include <savannah/systems/scene/scene.h>
// #include <savannah/systems/scene/sceneserializer.h>

// - ecs
// #include <savannah/systems/ecs/entity.h>
// #include <savannah/systems/ecs/components.h>

// - utils
// #include <savannah/utils/profiler.h>
// #include <savannah/utils/platform-utils.h> // modal windows

// entry point
// #include <savannah/main-loop.h> // need to be included in Client main file
#include <savannah/entry-point.h> // need to be included in Client main file

// Various Application data stored by Engine for Client App to use 
// #include <savannah/engine-data.h> // must be included last

/*
#define GL_SILENCE_DEPRECATION

// #include <savannah/external/glad/include/glad/glad.h>
#include <glfw/glfw3.h> // Will drag system OpenGL headers

// #define IMGUI_IMPL_OPENGL_LOADER_CUSTOM

#include <imgui/imgui.h>
#include <imgui/backends/imgui_impl_glfw.h>
#include <imgui/backends/imgui_impl_opengl3.h>
// #include <savannah/utils/dearimgui-backend.h>
*/

#endif
