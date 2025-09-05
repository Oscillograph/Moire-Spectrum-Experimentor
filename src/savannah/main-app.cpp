#include <savannah/main-app.h>

namespace Savannah
{
	App::App()
	{
	}
	
	App::~App()
	{
	}
	
	void App::Init()
	{
		PreSetup();
		SetupWindow();
		SetupResources();
		SetupGUI();
	}
	
	void App::PreSetup()
	{
	}
	
	void App::SetupWindow()
	{
		glfwSetErrorCallback(glfw_error_callback);
		if (!glfwInit())
			return;
		
		// Decide GL+GLSL versions
		// GL 3.0 + GLSL 130
		// const char* m_glsl_version = "#version 410";
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
		//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
		//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
		
		// Create window with graphics context
		// glfwWindowHint(GLFW_DECORATED, GLFW_FALSE);
		m_Window = glfwCreateWindow(m_WindowWidth, m_WindowHeight, m_WindowTitle.c_str(), nullptr, nullptr);
		if (m_Window == nullptr)
			return;
		glfwMakeContextCurrent(m_Window);
		gladLoadGL();
		glfwSwapInterval(1); // Enable vsync
	}
	
	void App::SetupResources()
	{
	}
	
	void App::SetupGUI()
	{
		// Setup Dear ImGui context
		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
        ImPlot::CreateContext();
		ImGuiIO& io = ImGui::GetIO(); // (void)io;
		io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
		io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
		io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;
		
		io.BackendFlags |= ImGuiBackendFlags_HasMouseCursors;
		io.BackendFlags |= ImGuiBackendFlags_HasSetMousePos;
		
		// Setup Dear ImGui style
		ImGui::StyleColorsDark();
		//ImGui::StyleColorsLight();
		
		// Setup Platform/Renderer backends
		ImGui_ImplOpenGL3_Init(m_glsl_version);
		ImGui_ImplGlfw_InitForOpenGL(m_Window, true);
		
		// load font
		if (font_config == nullptr){
			font_config = new ImFontConfig();
		}
		font_config->OversampleH = 1; //or 2 is the same
		font_config->OversampleV = 1;
		font_config->PixelSnapH = 1;
		io.Fonts->AddFontFromFileTTF("./data/fonts/tahoma.ttf", 18.0f, font_config, io.Fonts->GetGlyphRangesCyrillic());
		;
	}
	
	GLFWwindow* App::GetWindow()
	{
		return m_Window;
	}
	
	void App::ProcessInput()
	{
//		m_IdleFrames = 0;
		// Poll and handle events (inputs, window resize, etc.)
		// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
		// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
		// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
		// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
		UpdateWindowMinimizedStatus();
		glfwPollEvents();
	}
	
	void App::Logic()
	{
	}
	
	void App::Render()
	{
		if (!IsWindowMinimized())
		{
			// update current knowledge of the window size and scale
			glfwGetWindowSize(m_Window, &m_WindowCurrentWidth, &m_WindowCurrentHeight);
			m_WindowScale = {
				(float)m_WindowCurrentWidth / m_WindowWidth,
				(float)m_WindowCurrentHeight / m_WindowHeight 
			};
			
			GUIBegin();
			GUIContent();
			GUIEnd();
		}
	}
	
	void App::GUIBegin()
	{
		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
	}
	
	void App::GUIEnd()
	{
		// Rendering
		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(m_Window, &display_w, &display_h);
		// io.DisplaySize = ImVec2((float)display_w, (float)display_h);
		glViewport(0, 0, (GLsizei)display_w, (GLsizei)display_h);
		glClearColor(m_ClearColor.x * m_ClearColor.w, m_ClearColor.y * m_ClearColor.w, m_ClearColor.z * m_ClearColor.w, m_ClearColor.w);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		
		ImGuiIO& io = ImGui::GetIO();
		if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable){
			GLFWwindow* window_backup_context = glfwGetCurrentContext();
			ImGui::UpdatePlatformWindows();
			ImGui::RenderPlatformWindowsDefault();
			glfwMakeContextCurrent(window_backup_context);
		}
		
		// 
		glfwSwapBuffers(m_Window);
	}
	
	void App::GUIContent()
	{
		// Our state
		bool show_demo_window = true;
		bool show_another_window = false;
		m_ClearColor = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
		
		// 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
		if (show_demo_window)
			ImGui::ShowDemoWindow(&show_demo_window);
		// 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
		{
			ImGuiIO& io = ImGui::GetIO();
			static float f = 0.0f;
			static int counter = 0;
			
			ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.
			
			ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
			ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
			ImGui::Checkbox("Another Window", &show_another_window);
			
			ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
			ImGui::ColorEdit3("clear color", (float*)&m_ClearColor); // Edit 3 floats representing a color
			
			if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
				counter++;
			ImGui::SameLine();
			ImGui::Text("counter = %d", counter);
			
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
			ImGui::End();
		}
		
		// 3. Show another simple window.
		if (show_another_window)
		{
			ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
			ImGui::Text("Hello from another window!");
			if (ImGui::Button("Close Me"))
				show_another_window = false;
			ImGui::End();
		}
	}
	
	void App::ShutDown()
	{
		// Cleanup
		ImGui_ImplOpenGL3_Shutdown();
		ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
		ImGui::DestroyContext();
		
		glfwDestroyWindow(m_Window);
		glfwTerminate();
	}
	
	void App::SetFPS(float fps)
	{
		FPS = fps;
	}
	
	float App::GetFPS()
	{ 
		return FPS; 
	}
	
	int App::GetIdleFrames()
	{
		return m_IdleFrames;
	}
	
	void App::SetWindowTitle(const std::string& title)
	{
		m_WindowTitle = title;
	}
	
	void App::UpdateWindowMinimizedStatus()
	{
		int iconified = glfwGetWindowAttrib(GetWindow(), GLFW_ICONIFIED);
		if (iconified)
		{
			OnWindowMinimized();
		} else {
			OnWindowRestored();
		}
	}
	
	void App::OnWindowMinimized()
	{
		if (!m_WindowMinimized)
		{
			m_WindowMinimized = true;
		}
	}
	
	void App::OnWindowRestored()
	{
		if (m_WindowMinimized)
		{
			m_WindowMinimized = false;
		}
	}
	
	bool App::IsWindowMinimized()
	{
		return m_WindowMinimized;
	}
}
