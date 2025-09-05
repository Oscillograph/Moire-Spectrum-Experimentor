#ifndef SAVANNAH_MAINAPP_H
#define SAVANNAH_MAINAPP_H

#include <savannah/core.h>

// Dear ImGui: standalone example application for GLFW + OpenGL 3, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#include <stdio.h>

#include <external/glad/include/glad/glad.h>
#include <external/glfw/glfw3.h> // Will drag system OpenGL headers

#ifndef IMGUI_IMPL_OPENGL_LOADER_CUSTOM
#define IMGUI_IMPL_OPENGL_LOADER_CUSTOM
#endif

#include <external/imgui/imgui.h>
#include <external/imgui/imgui_internal.h> // for advanced selectable flags
#include <external/imgui/misc/cpp/imgui_stdlib.h> // for strings

#include <external/imgui/backends/imgui_impl_glfw.h>
#include <external/imgui/backends/imgui_impl_opengl3.h>

#include <external/imgui_filedialog/src/imgui_filedialog.h> // for file dialog
#include <external/implot/implot.h> // for plots

// #include <savannah/utils/dearimgui-backend.h>
#include <savannah/platforms/opengl/opengl_texture.h> // probably, it's better to make a separate renderer module

#define GL_SILENCE_DEPRECATION

static void glfw_error_callback(int error, const char* description)
{
	fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

namespace Savannah {
	class App {
	public:
		App();
		virtual ~App();
		virtual void Init();
		virtual void PreSetup();
		virtual void SetupWindow();
		virtual void SetupResources();
		virtual void SetupGUI();
		virtual GLFWwindow* GetWindow();
		virtual void ProcessInput();
		virtual void Logic();
		virtual void Render();
		virtual void GUIBegin();
		virtual void GUIEnd();
		virtual void GUIContent();
		virtual void ShutDown();
		virtual void SetFPS(float fps);
		virtual float GetFPS();
		virtual int GetIdleFrames();
		
		bool doExit = false;
		
	protected:
		ImFontConfig* font_config = nullptr;
		GLFWwindow* m_Window = nullptr;
		int m_WindowWidth = 1280;
		int m_WindowHeight = 720;
		int m_WindowCurrentWidth = 0;
		int m_WindowCurrentHeight = 0;
		ImVec2 m_WindowScale = {0.0f, 0.0f};
		ImVec4 m_ClearColor = {0.0f, 0.0f, 0.0f, 0.0f};
		std::string m_WindowTitle = "Превед, мир!";
		const char* m_glsl_version = "#version 410";
		bool m_WindowMinimized = false;
		float FPS = SAVANNAH_FPS60;
		int m_IdleFrames = 0;
		
		bool m_fileDialogOpen = false;
		ImFileDialogInfo m_fileDialogInfo;
		
	protected:
		void SetWindowTitle(const std::string& title);
		
	protected:
		void UpdateWindowMinimizedStatus();
		void OnWindowMinimized();
		void OnWindowRestored();
		bool IsWindowMinimized();
	};
	
	App* CreateApplication();
	
}

#endif
