#pragma once
#include "Observer.hpp"
#include "Utility.hpp"

#include "SPHSimulation.hpp"

struct GLFWwindow;

// TODO: Make singleton
class ImGuiViewer : public Observer {
public:
	ImGuiViewer() = default;
	~ImGuiViewer() override { Stop(); }

	void Start()
	{
		if (!m_Running.exchange(true))
			m_ImguiThread = std::thread([this]() { Loop(); });
	}
	void Stop ()
	{
		if (m_Running.exchange(false))
			m_ImguiThread.join();
	}

	/// Executes on the simulation thread
	void OnStartFrame() override;
	/// Executes on the simulation thread
	void OnEndFrame() override;
	/// Executes on the simulation thread
	void Attach(SPHSimulation* sim) override;

private:
	/// Executes on the UI thread
	void Loop();

	/// UI Setup
	bool Init();
	/// UI Cleanup
	void Deinit();

	/// Executes on the UI thread
	void DrawStatsWindow();
	/// Executes on the UI thread
	void DrawVisualizationWindow();

private:
	std::vector<Particle> m_Particles;
	SPHSimulation::grid_t m_Grid;

	GLFWwindow* window = nullptr;
	std::atomic<bool> m_Running = false;
	std::thread m_ImguiThread;
	std::mutex m_Mutex;

	bool m_Changed = false;

	float m_SimFPS = 0;
	stdc::time_point<stdclock> m_SimFrameStart;

	SPHSimulation::Params m_SimParams;

	enum ColoringParam
	{
		SUBDOMAIN, PRESSURE
	};
	ColoringParam m_ColoringParam = SUBDOMAIN;
};
