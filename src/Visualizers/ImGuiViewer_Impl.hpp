#ifndef DISABLE_UI
#include <format>

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
#include <imgui_internal.h>
#include <glad/glad.h>
#include <imgui_impl_opengl3.h>
#include <stdio.h>
#define GL_SILENCE_DEPRECATION
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <GLES2/gl2.h>
#endif

#include <GLFW/glfw3.h> // Will drag system OpenGL headers



#include "ImGuiWidgets.hpp"



static ImVec2 WorldToScreen(const Particle<2>::vec_t& p, const ImVec2 canvasPos, const ImVec2 canvasSize);
static Particle<2>::vec_t ScreenToWorld(const ImVec2&p, const ImVec2 canvasPos, const ImVec2 canvasSize);

static bool ComputeLinearRegression(const std::vector<coord<float, 2>>& points, float& outSlope, float& outIntercept);

template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::Attach(base::SPHSimulation<2, Particles>* sim)
{
	Observer<2, Particles>::Attach(sim);

	if (sim != nullptr)
	{
		m_SimParams = sim->GetParams();
		m_SimName = sim->GetName();
	}
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::OnEndFrame()
{
	if (this->m_Sim == nullptr)
		return;

	m_SimTimeCounter += stdclock::now() - m_SimFrameStart;
	m_SimFramesCounter++;
	if (m_SimTimeCounter > 500ms)
	{
		m_SimFPS = to<float>(m_SimFramesCounter) * 1'000'000'000.0 / m_SimTimeCounter.count();
		m_SimTrueFPS = to<float>(m_SimFramesCounter) * 1'000'000'000.0 / m_SimTrueTimeCounter.count();
		m_SimFramesCounter = 0;
		m_SimTimeCounter = 0ns;
		m_SimTrueTimeCounter = 0ns;
	}

	std::lock_guard<std::mutex> lock(m_Mutex);

	if (m_Changed)
	{
		this->m_Sim->SetParams(m_SimParams);
		this->m_Sim->ApplyCommand(m_Cmd);
		m_Changed = false;
	}

	// Copy, shouldn't be a bottleneck, and if it is then real time rendering isn't ideal anyway
	// TODO: there is no need for m_Particles to have the same memory layout of sim particles
	//       this could be turned in a copy like
	//       m_Sim->GetParticles().Populate(m_Particles)
	//       which fills the particles efficiently (for SoA resize and populate all positions first, then velocity...)
	
	if (m_RequestNewParticles)
	{
		m_RequestNewParticles = false;
		m_Particles = this->m_Sim->GetParticles();
	}
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::OnStartFrame()
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	if (this->m_Sim == nullptr)
		return;

	m_SimParams = this->m_Sim->GetParams();
	m_SimTime = this->m_Sim->GetTime();
	m_SimProfiling = this->m_Sim->GetProfiling();

	auto& objs = this->m_Sim->GetObjects();
	m_ObjectPositions.clear();
	m_ObjectTypes.clear();
	m_ObjectPositions.reserve(objs.size());
	m_ObjectTypes.reserve(objs.size());
	for (auto& obj : objs)
	{
		auto pos = obj->GetPosition();
		coord<float, 2> tl, br;
		tl.x = std::min(pos.first.x, pos.second.x);
		br.x = std::max(pos.first.x, pos.second.x);
		tl.y = std::max(pos.first.y, pos.second.y);
		br.y = std::min(pos.first.y, pos.second.y);
		this->m_ObjectPositions.push_back({ tl, br });
		this->m_ObjectTypes.push_back(obj->GetType());
	}

	auto now = stdclock::now();
	m_SimTrueTimeCounter += now - m_SimFrameStart;
	m_SimFrameStart = now;
}


template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::ProcessShortcuts()
{
	m_WantPasteProbes = ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_V);
	m_WantCopyProbes = ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_C);
}


template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::RenderStatsWindow()
{
#ifdef IMGUI_HAS_DOCK
	ImGui::SetNextWindowDockID(m_DockIDLeft, ImGuiCond_Appearing);
#endif
	ImGui::Begin("SPH Stats");

	ImGui::TextUnformatted(m_SimName.c_str());

	ImGui::Separator();
	{
		m_Changed |= ImGui::SliderFloat("Rest Density", &m_SimParams.RestDensity, 500, 1500, "%.0f");
		m_Changed |= ImGui::SliderFloat("Stiffness", &m_SimParams.Stiffness, 0, 1000, "%.0f");
		m_Changed |= ImGui::SliderFloat("Viscosity", &m_SimParams.Viscosity, 0, 5e-4, "%.5f");
		m_Changed |= ImGui::SliderFloat("Timestep", &m_SimParams.TimeStep, 0.0000f, 0.001f, "%.7f");
		m_Changed |= ImGui::SliderFloat("Smoothing Length", &m_SimParams.SmoothingLength, 0.0001f, 0.01f, "%.7f");
		m_Changed |= ImGui::SliderFloat("Final Time", &m_SimParams.FinalTime, 0.0f, 100.0f, "%.2f");
	}

	ImGui::Separator();
	{
		ImGui::Text("Particles: %ld", m_Particles.Size());
		ImGui::Text("UI / SPH (True) FPS: %.1f %.1f (%.1f)", ImGui::GetIO().Framerate, m_SimFPS, m_SimTrueFPS);

		{
			float progress = 0.0f;
			if (m_SimParams.FinalTime > 0.0f)
				progress = m_SimTime / m_SimParams.FinalTime;
			progress = ImClamp(progress, 0.0f, 1.0f);

			std::string label = std::format("{:.2f}/{:.2f}s", m_SimTime, m_SimParams.FinalTime);
			ImGui::ProgressBar(progress, ImVec2(-FLT_MIN, 0.0f), label.c_str());
			//ImGui::Text("Progress: %.2f / %.2f sec", m_SimTime, m_SimParams.FinalTime);
		}

		{
			std::vector<float> parts = { to<float>(m_SimProfiling.Neighbors.count()), to<float>(m_SimProfiling.Initialize.count()), to<float>(m_SimProfiling.IterativePressure.count()) };
			int hovering = DrawStackedProgressBar(parts);
			if (hovering != -1)
			{
				static std::string labels[] = { "Neighbors", "Initialize", "Iterative Pressure" };

				float sum = 0.0f;
				for (float v : parts)
					sum += ImMax(v, 0.0f);

				ImGui::BeginTooltip();
				std::string label = std::format("{}: {:.1f}%%", labels[hovering], parts[hovering] / sum * 100);
				ImGui::TextUnformatted(label.c_str());
				ImGui::EndTooltip();
			}
		}
	}
	
	ImGui::Separator();
	if (!m_Particles.Empty())
	{
		ImGui::Text("Position: %.3f %.3f", m_Particles.PositionX(0), m_Particles.PositionY(0));
		ImGui::Text("Velocity: %.3f %.3f", m_Particles.VelocityX(0), m_Particles.VelocityY(0));
	}
	if (m_Cmd.Type != Command<2>::NONE)
		ImGui::Text("Command: %.3f %.3f %.3f %.3f", m_Cmd.Position.x, m_Cmd.Position.y, m_Cmd.Radius, m_Cmd.Strength);
	else
		ImGui::Text("Command: NONE");

	ImGui::Separator();
	{
		static const char* items[] = { "Pressure", "Velocity", "Density"};
		static int coloring_param_int = 1;
		ImGui::Combo("Coloring Parameter", &coloring_param_int, items, IM_ARRAYSIZE(items));
		m_ColoringParam = (ColoringParam)coloring_param_int;

		if (m_ColoringParam == PRESSURE)
		{
			ImGui::SliderFloat("Max UI Pressure", &m_MaxPressure, 0, 20000, "%.0f");
			DrawLegendScale(ImVec2(200.0f, ImGui::GetTextLineHeightWithSpacing()), 0, m_MaxPressure);
		}
		else if (m_ColoringParam == VELOCITY)
		{
			ImGui::SliderFloat("Max UI Velocity", &m_MaxVelocity, 0, 5, "%.5f");
			DrawLegendScale(ImVec2(200.0f, ImGui::GetTextLineHeightWithSpacing()), 0, m_MaxVelocity);
		}
		else if (m_ColoringParam == DENSITY)
		{
			ImGui::SliderFloat("Max UI Density", &m_MaxDensity, 0, 10000, "%.5f");
			DrawLegendScale(ImVec2(200.0f, ImGui::GetTextLineHeightWithSpacing()), 0, m_MaxDensity);
		}
	}

	ImGui::Separator();
	{
		if (ImGui::Button("Probes"))
			m_ShowProbes = true;
	}


	ImGui::End();
}


template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::RenderVisualizationWindow()
{
#ifdef IMGUI_HAS_DOCK
	ImGui::SetNextWindowDockID(m_DockIDCenter, ImGuiCond_Appearing);
#endif
	ImGui::Begin("Fluid");

	ImVec2 canvasPos  = ImGui::GetCursorScreenPos();
	ImVec2 canvasSize = ImGui::GetContentRegionAvail();


	// Click
	ImVec2 mousePos = ImGui::GetMousePos();
	vec_t worldPos  = ScreenToWorld(mousePos, canvasPos, canvasSize);

	bool hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
	bool clicked = ImGui::IsMouseDown(ImGuiMouseButton_Left) && worldPos.x > 0 && worldPos.y > 0 && worldPos.x < 1.0 && worldPos.y < 1.0;
	if (hovered && clicked)
	{
		m_Cmd.Type = Command<2>::PRESSURE;
		m_Cmd.Position = worldPos;
		m_Cmd.Radius = 0.1;
		m_Cmd.Strength = 20000;
		m_Changed = true;
	}
	else if (m_Cmd.Type != Command<2>::NONE)
	{
		m_Cmd.Type = Command<2>::NONE;
		m_Changed = true;
	}

	// Draw
	ImDrawList* drawList = ImGui::GetWindowDrawList();

	// Background
	drawList->AddRectFilled(
		canvasPos,
		ImVec2(canvasPos.x + canvasSize.x, canvasPos.y + canvasSize.y),
		IM_COL32(20, 20, 20, 255)
	);

	DrawParticles(drawList);
	DrawObjects(drawList);
	DrawProbes(drawList);

	ImGui::End();
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::DrawParticles(ImDrawList* drawList)
{
	ImVec2 canvasPos = ImGui::GetCursorScreenPos();
	ImVec2 canvasSize = ImGui::GetContentRegionAvail();

	for (size_t i = 0; i < m_Particles.Size(); i++) {
		const Particle<2>& p = m_Particles.GetParticle(i);

		ImVec2 pos = WorldToScreen(p.Position, canvasPos, canvasSize);
		static constexpr float r = 2.0f;
		ImU32 color = ImColor::HSV(1.0f, 1.0f, 1.0f);

		if (p.Type == SOLID)
		{
			color = ImColor(255, 255, 255);
		}
		else if (m_ColoringParam == PRESSURE)
		{
			float t = std::clamp(p.Pressure / m_MaxPressure, 0.0f, 1.0f);
			color = ImColor::HSV(0.66f - 0.66f * t, 1.0f, 1.0f); // blue -> red
		}
		else if (m_ColoringParam == DENSITY)
		{
			float t = std::clamp(p.Density / m_MaxDensity, 0.0f, 1.0f);
			color = ImColor::HSV(0.66f - 0.66f * t, 1.0f, 1.0f); // blue -> red
		}
		else if (m_ColoringParam == VELOCITY)
		{
			float t = Norm(p.Velocity);
			t = std::clamp(t / m_MaxVelocity, 0.0f, 1.0f);
			color = ImColor::HSV(0.66f - 0.66f * t, 1.0f, 1.0f); // blue -> red
		}

		drawList->AddCircleFilled(pos, r, color);
	}
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::DrawObjects(ImDrawList* drawList)
{
	ImVec2 canvasPos = ImGui::GetCursorScreenPos();
	ImVec2 canvasSize = ImGui::GetContentRegionAvail();

	for (int i = 0; i < m_ObjectPositions.size(); i++)
	{
		auto pos = m_ObjectPositions[i];
		auto type = m_ObjectTypes[i];

		ImU32 color;

		if (type == ObjectType::UNKNOWN)
			color = ImColor(216, 0, 255);
		else if (type == ObjectType::SOURCE)
			color = ImColor(83, 154, 255);
		else if (type == ObjectType::SINK)
			color = ImColor(255, 50, 50);

		ImVec2 tl = WorldToScreen(pos.first, canvasPos, canvasSize), br = WorldToScreen(pos.second, canvasPos, canvasSize);

		drawList->AddRect(tl, br, color);
		drawList->AddLine(tl, br, color);
	}
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::DrawProbes(ImDrawList* drawList)
{
	ImVec2 canvasPos = ImGui::GetCursorScreenPos();
	ImVec2 canvasSize = ImGui::GetContentRegionAvail();

	for (int i = 0; i < m_Probes.size(); i++)
	{
		Probe& p = m_Probes[i];

		ImVec2 tl = WorldToScreen(p.TL, canvasPos, canvasSize), br = WorldToScreen(p.BR, canvasPos, canvasSize);
		ImU32 color = ImColor(255, 255, 255);

		drawList->AddRect(tl, br, color);
		drawList->AddLine(tl, br, color);
	}
}


template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::RenderProbeWindow()
{
	ImGuiContext& g = *GImGui;
	const ImGuiStyle& style = g.Style;
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 center = ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f);

	static stdclock::time_point error_time;
	static std::string error_msg = "";
	static bool error_show = false;

	if (!m_ShowProbes)
		return;

#ifdef IMGUI_HAS_DOCK
	ImGui::SetNextWindowDockID(m_DockIDBottom, ImGuiCond_Appearing);
#endif
	bool visible = ImGui::Begin("Probes", &m_ShowProbes);
	if (!visible)
	{
		ImGui::End();
		return;
	}

	// Add Probe
	{
		static bool show_edit_window = false;
		static float p1[2] = { 0.0f, 0.0f }, p2[2] = { 0.0f, 0.0f };
		static constexpr size_t SIZE = 255;
		static char buf[SIZE] = "\0";

		if (m_EditingProbe)
			ImGui::BeginDisabled();

		if (ImGui::Button("Add"))
		{
			show_edit_window = true;
			p1[0] = 0; p1[1] = 0; p2[0] = 0; p2[1] = 0;
			buf[0] = '\0';
		}

		if (m_EditingProbe)
			ImGui::EndDisabled();

		// use a static variable to show the window m_EditingProbe can be set from other loc
		if (show_edit_window)
		{
			m_EditingProbe = true;
			if (BeginEditProbeWindow(center, &m_EditingProbe, p1, p2, buf, SIZE, m_WantCopyProbes, m_WantPasteProbes))
			{
				Probe p;
				p.Name = buf;
				if (p.Name.empty())
					p.Name = std::format("Probe {}", m_Probes.size() + 1);
				p.Selected = false;
				p.TL = { p1[0], p1[1] };
				p.BR = { p2[0], p2[1] };
				m_Probes.push_back(p);
				m_EditingProbe = false;
				show_edit_window = false;
			}
			if (m_EditingProbe == false)
				show_edit_window = false;
		}
	}
			
	ImGui::SameLine();
	if (ImGui::Button("Remove All"))
		m_Probes.clear();
			
	ImGui::SameLine();

	// Select All
	{
		bool are_all_selected = !m_Probes.empty();
		for (auto& p : m_Probes)
			if (!p.Selected)
				are_all_selected = false;

		if (are_all_selected)
		{
			if (ImGui::Button("Deselect All"))
			{
				for (auto& p : m_Probes)
					p.Selected = false;
			}
		}
		else
		{
			if (ImGui::Button("Select All"))
			{
				for (auto& p : m_Probes)
					p.Selected = true;
			}
		}
	}

	// Copy/Paste
	{
		if (ImGui::Button("Copy") || (!m_EditingProbe && m_WantCopyProbes))
		{
			std::stringstream ss;
			
			for (const auto& p : m_Probes)
			{
				ss << p.ToString();
				ss << ';';
			}

			ImGui::SetClipboardText(ss.str().c_str());
		}
		ImGui::SameLine();

		if (ImGui::Button("Paste") || (!m_EditingProbe && m_WantPasteProbes))
		{
			const char* clip = ImGui::GetClipboardText();
			if (clip)
			{
				std::string str = clip;
				std::optional<std::vector<Probe>> op = Probe::ParseMultiple(str);
				if (op.has_value())
				{
					std::vector<Probe> ps = op.value();
					m_Probes.insert(m_Probes.end(), ps.begin(), ps.end());
				}
				else
				{
					error_msg = "Invalid text pasted, use format {{name}:{tl.x}:{tl.y}:{br.x}:{br.y}};{...}";
					error_time = stdclock::now();
					error_show = true;
				}
			}
		}

		ImGui::SameLine();

		static const ImGuiComboFlags flags = ImGuiComboFlags_None;
		static const char* items[] = { "Box Small", "Box Big" };
		static const Probe probes_list[] = { Probe::ParseOne("Probe Small:0.5:0.1:0.51:0.045").value(), Probe::ParseOne("Probe Big:0.5:0.5:0.51:0.045").value() };
		if (ImGui::BeginCombo("##CommonProbes", "Common Probes", flags))
		{
			for (int n = 0; n < IM_ARRAYSIZE(items); n++)
				if (ImGui::Selectable(items[n], false))
					m_Probes.push_back(probes_list[n]);
			ImGui::EndCombo();
		}

	}
			
	{
		if (ImGui::Button("Compute"))
			m_ShowProbesData = true;

		if (m_ShowProbesData)
			RenderProbeDataWindow();
	}

	if (error_show)
	{
		ImGui::TextUnformatted(error_msg.c_str());

		if (stdclock::now() - error_time > 2s)
			error_show = false;
	}

	if (ImGui::BeginChild("ProbesList", ImVec2(0, 0), true))
	{
		for (size_t i = 0; i < m_Probes.size(); ++i)
		{
			if (RenderProbe(i))
			{
				i--;
				break;
			}
		}
	}
	ImGui::EndChild();

	ImGui::End();
}
template <ParticleSet<2> Particles>
bool ImGuiViewer<Particles>::RenderProbe(size_t i)
{
	ImGuiContext& g = *GImGui;
	const ImGuiStyle& style = g.Style;
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 center = ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f);
	ImVec2 padding = ImGui::GetStyle().FramePadding;

	static int edit = -1;

	ImGui::PushID(to<int>(i)); // ensure unique IDs

	auto& probe = m_Probes[i];

	// Checkbox
	ImGui::Checkbox(probe.Name.c_str(), &probe.Selected);

	// Right-aligned "Remove" button
	float removeWidth = ImGui::CalcTextSize("Remove").x + padding.x * 2;
	float editWidth = ImGui::CalcTextSize("Edit").x + padding.x * 2;
	float avail = ImGui::GetContentRegionAvail().x;

	ImGui::SameLine(avail - editWidth - removeWidth - padding.x);

	bool should_disable = m_EditingProbe;

	// Edit
	{
		static float p1[2] = { 0.0f, 0.0f }, p2[2] = { 0.0f, 0.0f };
		static constexpr size_t SIZE = 255;
		static char buf[SIZE] = "\0";

		if (should_disable)
			ImGui::BeginDisabled();

		if (ImGui::Button("Edit"))
		{
			edit = i;
			m_EditingProbe = true;
			p1[0] = probe.TL.x; p1[1] = probe.TL.y; p2[0] = probe.BR.x; p2[1] = probe.BR.y;
			strncpy(buf, probe.Name.substr(0, SIZE - 1).data(), SIZE - 1);
			buf[SIZE - 1] = '\0';
		}

		if (should_disable)
			ImGui::EndDisabled();

		if (edit == i)
		{
			if (BeginEditProbeWindow(center, &m_EditingProbe, p1, p2, buf, SIZE, m_WantCopyProbes, m_WantPasteProbes))
			{
				probe.Name = buf;
				if (probe.Name.empty())
					probe.Name = std::format("Probe {}", m_Probes.size() + 1);
				probe.TL = { p1[0], p1[1] };
				probe.BR = { p2[0], p2[1] };
				edit = -1;
				m_EditingProbe = false;
			}
			if (m_EditingProbe == false)
				edit = -1;
		}
	}

	ImGui::SameLine(avail - removeWidth);

	if (should_disable)
		ImGui::BeginDisabled();

	if (ImGui::Button("Remove"))
	{
		m_Probes.erase(m_Probes.begin() + i);

		if (edit != -1)
		{
			edit = -1;
			m_EditingProbe = false;
		}

		ImGui::PopID();
		return true;
	}

	if (should_disable)
		ImGui::EndDisabled();

	ImGui::PopID();
	return false;
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::RenderProbeDataWindow()
{
	ImGuiContext& g = *GImGui;
	const ImGuiStyle& style = g.Style;
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 center = ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f);

	ImGui::SetNextWindowPos(center, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(800, 400), ImGuiCond_FirstUseEver);

	bool visible = ImGui::Begin("Data", &m_ShowProbesData);
	if (!visible) {
		ImGui::End();
		return;
	}

	enum class PlotType : int {
		PRESSURE_DEPTH = 0
	};

	static const char* items[] = { "Pressure/Depth" };
	static PlotType plot_type = PlotType::PRESSURE_DEPTH;
	ImGui::TextUnformatted("Plot Type: "); ImGui::SameLine();
	ImGui::Combo("##PlotType", (int*)&plot_type, items, IM_ARRAYSIZE(items));

	ImGui::SameLine();

	/// points sampled from a probe should be contiguous
	static std::vector<coord<float, 2>> points;
	/// number of points sampled from the i-th probe, and it's ID
	static std::vector<std::pair<int, size_t>> counts;

	static bool realtime = false;
	if (ImGui::Button("Generate") || realtime)
	{
		// Sample data
		points.clear();

		if (plot_type == PlotType::PRESSURE_DEPTH)
			SampleForPressureDepth(points, counts);
	}

	ImGui::SameLine();
	ImGui::Checkbox("Realtime", &realtime);

	ImVec2 canvasPos = ImGui::GetCursorScreenPos();
	ImVec2 canvasSize = ImGui::GetContentRegionAvail();
	canvasSize.y = canvasSize.y - ImGui::GetFrameHeightWithSpacing();

	GraphRange range;

	if (plot_type == PlotType::PRESSURE_DEPTH)
		range = DrawPressureDepthGraph(points, counts, canvasSize, std::get<AdditionalPressureDepthSampleData>(m_AdditionalSampleData));

	ImGui::Text("Items: %ld", points.size());
	
	ImVec2 mousePos = ImGui::GetMousePos();
	vec_t worldPos = ScreenToWorld(mousePos, canvasPos, canvasSize);
	if (worldPos.x > 0 && worldPos.y > 0 && worldPos.x < 1.0 && worldPos.y < 1.0)
	{
		ImGui::SameLine();
		ImGui::Text("| Mouse Position: %.3f:%.3f", worldPos.x * (range.x.max - range.x.min) + range.x.min, worldPos.y * (range.y.max - range.y.min) + range.y.min);
	}

	ImGui::End();
}


template <ParticleSet<2> Particles>
GraphRange ImGuiViewer<Particles>::DrawPressureDepthGraph(const std::vector<coord<float, 2>>& points, const std::vector<std::pair<int, size_t>>& counts, ImVec2 size, const AdditionalPressureDepthSampleData& data)
{
	auto graph_range = DrawGraph(points, size, [&](ImDrawList* drawList, ImVec2 min, ImVec2 max, ImVec2 range, ImVec2 canvasPos, ImVec2 canvasSize)
		{
			for (const auto& p : points)
				drawList->AddCircleFilled(ToScreenSpace(p, min, range, canvasPos, canvasSize), 3.0f, IM_COL32(255, 255, 255, 255));
			
			if (data.valid)
			{
				coord<float, 2> p1{ min.x, data.m * min.x + data.b };
				coord<float, 2> p2{ max.x, data.m * max.x + data.b };

				ImVec2 sp1 = ToScreenSpace(p1, min, range, canvasPos, canvasSize);
				ImVec2 sp2 = ToScreenSpace(p2, min, range, canvasPos, canvasSize);

				drawList->AddLine(sp1, sp2, IM_COL32(255, 0, 0, 255), 2.0f);
			}

			std::string text = std::format("Slope = {}; Intercept = {}", data.m, data.b);
			drawList->AddText(ImVec2{ canvasPos.x + 10, canvasPos.y + 10 }, IM_COL32(255, 255, 255, 255), text.c_str());
		});

	return graph_range;
}


template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::BeginFullscreenDockspace()
{
#ifdef IMGUI_HAS_DOCK
	const ImGuiViewport* viewport = ImGui::GetMainViewport();
	ImGui::DockSpaceOverViewport(m_DockspaceID, viewport, ImGuiDockNodeFlags_PassthruCentralNode);
#endif
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::BuildInitialLayout()
{
#ifdef IMGUI_HAS_DOCK
	static bool first_time = true;
	if (!first_time) return;
	first_time = false;

	m_DockspaceID = ImGui::GetID("FullscreenDockspace");
	ImGui::DockBuilderRemoveNode(m_DockspaceID);
	ImGui::DockBuilderAddNode(m_DockspaceID, ImGuiDockNodeFlags_DockSpace);
	ImGui::DockBuilderSetNodeSize(m_DockspaceID, ImVec2(800, 600));

	m_DockIDLeft = ImGui::DockBuilderSplitNode(m_DockspaceID, ImGuiDir_Left, 0.40f, nullptr, &m_DockIDCenter);
	m_DockIDBottom = ImGui::DockBuilderSplitNode(m_DockIDLeft, ImGuiDir_Down, 0.30f, nullptr, &m_DockIDLeft);

	ImGui::DockBuilderFinish(m_DockspaceID);
#endif
}


template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::ForAllParticlesInProbes(const std::function<void(int, int, float, float)>& fn)
{
	for (const Probe& probe : m_Probes)
	{
		if (!probe.Selected)
			continue;

		for (size_t i = 0; i < m_Particles.Size(); i++)
		{
			if (m_Particles.Type(i) != ParticleType::FLUID)
				continue;

			coord<float, 2> pos = m_Particles.Position(i);

			if (pos.x >= probe.TL.x && pos.x <= probe.BR.x && pos.y >= probe.BR.y && pos.y <= probe.TL.y)
				fn(probe.ID, i, pos.x, pos.y);
		}
	}
}
template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::SampleForPressureDepth(std::vector<coord<float, 2>>& out_points, std::vector<std::pair<int, size_t>>& out_counts)
{
	struct Data
	{
		float y, pressure;
	};

	std::vector<Data> data;
	int last_probe_id = Probe::INVALID_ID;
	size_t count = 0;

	ForAllParticlesInProbes([&](int probe_id, int particle_index, float x, float y)
		{
			data.push_back({ y, m_Particles.Pressure(particle_index) });
			
			if (probe_id != last_probe_id)
			{
				if (last_probe_id != Probe::INVALID_ID)
					out_counts.push_back({ last_probe_id, count });
				count = 1;
				last_probe_id = probe_id;
			}
			else
				count++;
		});

	if (last_probe_id != Probe::INVALID_ID)
		out_counts.push_back({ last_probe_id, count });

	if (data.empty())
		return;

	float maxy = data[0].y;
	for (const Data& part : data)
		maxy = std::max(maxy, part.y);

	out_points.reserve(data.size());
	for (const Data& part : data)
		out_points.push_back({ -part.y + maxy, part.pressure });

	AdditionalPressureDepthSampleData add_data;
	add_data.valid = ComputeLinearRegression(out_points, add_data.m, add_data.b);
	m_AdditionalSampleData = add_data;
}




static ImVec2 WorldToScreen(const Particle<2>::vec_t& p, const ImVec2 canvasPos, const ImVec2 canvasSize) {
	return ImVec2(
		canvasPos.x + p.x * canvasSize.x,
		canvasPos.y + (1.0f - p.y) * canvasSize.y
	);
}
static Particle<2>::vec_t ScreenToWorld(const ImVec2& mouse, const ImVec2 canvasPos, const ImVec2 canvasSize) {
	return Particle<2>::vec_t{
			   to<float>(mouse.x - canvasPos.x) / canvasSize.x,
		1.0f - to<float>(mouse.y - canvasPos.y) / canvasSize.y
	};
}


static void glfw_error_callback(int error, const char* description)
{
	fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::Loop()
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

		static bool show_demo_window = false;
		static ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

		if (show_demo_window)
			ImGui::ShowDemoWindow(&show_demo_window);

		BuildInitialLayout();
		BeginFullscreenDockspace();

		{
			std::lock_guard<std::mutex> lock(m_Mutex);

			ProcessShortcuts();
			RenderStatsWindow();
			RenderVisualizationWindow();
			RenderProbeWindow();

			m_RequestNewParticles = true;
		}

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


template <ParticleSet<2> Particles>
void ImGuiViewer<Particles>::Deinit()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();
}
template <ParticleSet<2> Particles>
bool ImGuiViewer<Particles>::Init()
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
	window = glfwCreateWindow((int)(1280 * main_scale), (int)(800 * main_scale), "Cool 2D SPH visualization", nullptr, nullptr);
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
	#ifdef IMGUI_HAS_DOCK
	io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
	#endif

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

	io.ConfigWindowsMoveFromTitleBarOnly = true;

	return true;
}
#endif



static bool ComputeLinearRegression(const std::vector<coord<float, 2>>& points, float& outSlope, float& outIntercept)
{
	size_t n = points.size();
	if (n < 2)
		return false;

	double sumX = 0.0;
	double sumY = 0.0;
	double sumXY = 0.0;
	double sumXX = 0.0;

	for (const auto& p : points)
	{
		sumX += p.x;
		sumY += p.y;
		sumXY += p.x * p.y;
		sumXX += p.x * p.x;
	}

	double denom = (n * sumXX - sumX * sumX);

	// Prevent division by zero (vertical line or identical Xs)
	if (std::abs(denom) < 1e-8)
		return false;

	outSlope = static_cast<float>((n * sumXY - sumX * sumY) / denom);
	outIntercept = static_cast<float>((sumY - outSlope * sumX) / n);

	return true;
}