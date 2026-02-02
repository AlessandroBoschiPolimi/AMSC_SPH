#pragma once
#include "Observer.hpp"
#include "Utility.hpp"

#include "SPHSimulation.hpp"

#include "Command.hpp"

struct GLFWwindow;

// TODO: Make singleton
class ImGuiViewer : public Observer<2> {
public:
	using vec_t = Particle<2>::vec_t;

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
	void Attach(SPHSimulation<2>* sim) override;

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
	std::vector<Particle<2>> m_Particles;
	SPHSimulation<2>::grid_t m_Grid;

	GLFWwindow* window = nullptr;
	std::atomic<bool> m_Running = false;
	std::thread m_ImguiThread;
	std::mutex m_Mutex;

	bool m_Changed = false;

	float m_SimFPS = 0;
	stdc::time_point<stdclock> m_SimFrameStart;
	float m_MaxVelocity = 10, m_MaxPressure = 50000;
	Command<2> m_Cmd;

	SPHSimulation<2>::Params m_SimParams;

	enum ColoringParam
	{
		SUBDOMAIN, PRESSURE, VELOCITY
	};
	ColoringParam m_ColoringParam = VELOCITY;
};
