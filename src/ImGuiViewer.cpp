#include "ImGuiViewer.hpp"
#include "SPHSimulation.hpp"

// Dear ImGui: standalone example application for GLFW + OpenGL 3, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)

// Learn about Dear ImGui:
// - FAQ                  https://dearimgui.com/faq
// - Getting Started      https://dearimgui.com/getting-started
// - Documentation        https://dearimgui.com/docs (same as your local docs/ folder).
// - Introduction, links and more at the top of imgui.cpp

#define IMGUI_IMPL_OPENGL_LOADER_GLAD
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <glad/glad.h>
#include <imgui_impl_opengl3.h>
#include <stdio.h>
#define GL_SILENCE_DEPRECATION
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <GLES2/gl2.h>
#endif

#include <GLFW/glfw3.h> // Will drag system OpenGL headers


// [Win32] Our example includes a copy of glfw3.lib pre-compiled with VS2010 to maximize ease of testing and compatibility with old VS compilers.
// To link with VS2010-era libraries, VS2015+ requires linking with legacy_stdio_definitions.lib, which we do using this pragma.
// Your own project should not be affected, as you are likely to link with a newer binary of GLFW that is adequate for your version of Visual Studio.
#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif



void ImGuiViewer::Attach(SPHSimulation* sim)
{
	Observer::Attach(sim);

	m_SimParams = sim->GetParams();
}
void ImGuiViewer::OnEndFrame()
{
	auto simFramEnd = stdclock::now();
	auto dt = simFramEnd - m_SimFrameStart;
	m_SimFPS = 1'000'000'000.0 / dt.count();

	// TODO: ASSERT m_Sim != nullptr
	std::lock_guard<std::mutex> lock(m_Mutex);

	if (m_Changed)
	{
		m_Sim->SetParams(m_SimParams);
	}

	// Copy, shouldn't be a bottleneck, and if it is then real time rendering isn't ideal anyway
	// TODO: triple buffering
	m_Particles = m_Sim->GetParticles();
	m_Grid = m_Sim->GetGrid();
}
void ImGuiViewer::OnStartFrame()
{
	// TODO: ASSERT m_Sim != nullptr
	std::lock_guard<std::mutex> lock(m_Mutex);

	m_SimParams = m_Sim->GetParams();

	m_SimFrameStart = stdclock::now();
}



void ImGuiViewer::DrawStatsWindow() {
	ImGui::Begin("SPH Stats");

	{
		std::lock_guard<std::mutex> lock(m_Mutex);

		m_Changed |= ImGui::SliderFloat("Rest Density", &m_SimParams.RestDensity, 500      , 1500     , "%.5f");
		m_Changed |= ImGui::SliderFloat("Stiffness"   , &m_SimParams.Stiffness  ,   0      , 5000     , "%.5f");
		m_Changed |= ImGui::SliderFloat("Viscosity"   , &m_SimParams.Viscosity  ,   0      ,    1     , "%.5f");
		m_Changed |= ImGui::SliderFloat("dt"          , &m_SimParams.TimeStep   ,   0.0001f,    0.005f, "%.5f");
		m_Changed |= ImGui::SliderFloat("Smoothing Length", &m_SimParams.SmoothingLength,   0.0001f,    0.5f, "%.5f");

		ImGui::Text("Particles: %d", m_Particles.size());
		ImGui::Text("UI  FPS: %.1f", ImGui::GetIO().Framerate);
		ImGui::Text("SPH FPS: %.1f", m_SimFPS);

		const char* items[] = { "Subdomain", "Pressure" };
		static int coloring_param_int = 0;
		ImGui::Combo("Coloring Parameter", &coloring_param_int, items, IM_ARRAYSIZE(items));
		m_ColoringParam = (ColoringParam)coloring_param_int;
	}

	ImGui::End();
}
void ImGuiViewer::DrawVisualizationWindow() {
	ImGui::Begin("Fluid");

	ImVec2 canvasPos = ImGui::GetCursorScreenPos();
	ImVec2 canvasSize = ImGui::GetContentRegionAvail();

	ImDrawList* drawList = ImGui::GetWindowDrawList();

	// Background
	drawList->AddRectFilled(
		canvasPos,
		ImVec2(canvasPos.x + canvasSize.x, canvasPos.y + canvasSize.y),
		IM_COL32(20, 20, 20, 255)
	);

	// World -> screen transform
	auto worldToScreen = [&](const coord<float, 2>& p) {
		return ImVec2(
			canvasPos.x + p.x * canvasSize.x,
			canvasPos.y + (1.0f - p.y) * canvasSize.y
		);
		};

	// Draw particles
	for (const auto& p : m_Particles) {
		ImVec2 pos = worldToScreen(p.Position);
		float r = 2.0f;
		ImU32 color = ImColor::HSV(1.0f, 1.0f, 1.0f);
		if (p.Type == SOLID){
			color = ImColor::HSV(
				0.0f, // blue -> red
				0.0f,
				1.0f
			);
			drawList->AddCircleFilled(pos, r, color);
			continue;
		}
		if (m_ColoringParam == PRESSURE)
		{
			float d = p.Density;
			float t = (d - 0.1) / (1.0 - 0.1);
			t = std::clamp(t, 0.0f, 1.0f);

			color = ImColor::HSV(
				0.6f - 0.6f * t, // blue -> red
				1.0f,
				1.0f
			);
		}
		else if (m_ColoringParam == SUBDOMAIN)
		{
			SPHSimulation::cell_pos_t cell_pos = SPHSimulation::GetCellPosition(p.Position, m_SimParams.SmoothingLength);
			color = ImColor::HSV(
				1.0f - cell_pos.x * 0.1,
				1.0f - cell_pos.y * 0.1,
				1.0f
			);
		}

		drawList->AddCircleFilled(pos, r, color);
	}

	ImGui::End();
}


static void glfw_error_callback(int error, const char* description)
{
	fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

void ImGuiViewer::Loop()
{
	if (!Init())
	{
		m_Running = false;
		return;
	}

	while (m_Running && !glfwWindowShouldClose(window))
	{
		// Poll and handle events (inputs, window resize, etc.)
		// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
		// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
		// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
		// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
		glfwPollEvents();
		if (glfwGetWindowAttrib(window, GLFW_ICONIFIED) != 0)
		{
			ImGui_ImplGlfw_Sleep(10);
			continue;
		}

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		static bool show_demo_window = true;
		static ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

		if (show_demo_window)
			ImGui::ShowDemoWindow(&show_demo_window);

		DrawStatsWindow();
		DrawVisualizationWindow();

		// Rendering
		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	Deinit();
	m_Running = false;
}


void ImGuiViewer::Deinit()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();
}
bool ImGuiViewer::Init()
{
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit())
		return false;

	// Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
	// GL ES 2.0 + GLSL 100 (WebGL 1.0)
	const char* glsl_version = "#version 100";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(IMGUI_IMPL_OPENGL_ES3)
	// GL ES 3.0 + GLSL 300 es (WebGL 2.0)
	const char* glsl_version = "#version 300 es";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
	// GL 3.2 + GLSL 150
	const char* glsl_version = "#version 150";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
	// GL 3.0 + GLSL 130
	const char* glsl_version = "#version 130";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

	// Create window with graphics context
	float main_scale = ImGui_ImplGlfw_GetContentScaleForMonitor(glfwGetPrimaryMonitor()); // Valid on GLFW 3.3+ only
	window = glfwCreateWindow((int)(1280 * main_scale), (int)(800 * main_scale), "Dear ImGui GLFW+OpenGL3 example", nullptr, nullptr);
	if (window == nullptr)
		return false;
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync

	gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsLight();

	// Setup scaling
	ImGuiStyle& style = ImGui::GetStyle();
	style.ScaleAllSizes(main_scale);        // Bake a fixed style scale. (until we have a solution for dynamic style scaling, changing this requires resetting Style + calling this again)
	style.FontScaleDpi = main_scale;        // Set initial font scale. (using io.ConfigDpiScaleFonts=true makes this unnecessary. We leave both here for documentation purpose)

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	// Load Fonts
	// - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
	// - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
	// - If the file cannot be loaded, the function will return a nullptr. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
	// - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
	// - Read 'docs/FONTS.md' for more instructions and details. If you like the default font but want it to scale better, consider using the 'ProggyVector' from the same author!
	// - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
	// - Our Emscripten build process allows embedding fonts to be accessible at runtime from the "fonts/" folder. See Makefile.emscripten for details.
	//style.FontSizeBase = 20.0f;
	//io.Fonts->AddFontDefault();
	//io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf");
	//io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf");
	//io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf");
	//io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf");
	//ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf");
	//IM_ASSERT(font != nullptr);

	return true;
}
