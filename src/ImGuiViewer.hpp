#pragma once
#include "Observer.hpp"
#include "Utility.hpp"

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

	void OnStartFrame() override;
	void OnEndFrame() override;
	void Attach(SPHSimulation* sim) override;

private:
	void Loop();

	void Deinit();
	bool Init();

	void DrawStatsWindow();
	void DrawVisualizationWindow();

private:
	std::vector<Particle> m_Particles;

	GLFWwindow* window = nullptr;
	std::atomic<bool> m_Running = false;
	std::thread m_ImguiThread;
	std::mutex m_Mutex;

	bool m_Changed = false;

	float m_SimFPS = 0;
	stdc::time_point<stdclock> m_SimFrameStart;

	float m_RestDensity = 0;
	float m_Stiffness   = 0;
	float m_Viscosity   = 0;
	float m_TimeStep    = 0;
	float m_Increment   = 0;
};