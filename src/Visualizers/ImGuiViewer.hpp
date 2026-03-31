#pragma once
#include "Observer.hpp"
#include "Utility.hpp"

#include "Base/SPHSimulation.hpp"

#include "Command.hpp"

struct GLFWwindow;

// TODO: Make singleton
template <ParticleSet<2> Particles>
class ImGuiViewer : public Observer<2, Particles> {
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
		if (m_ImguiThread.joinable())
		{
			m_Running = false;
			m_ImguiThread.join();
		}
	}

	/// Executes on the simulation thread
	void OnStartFrame() override;
	/// Executes on the simulation thread
	void OnEndFrame() override;
	/// Executes on the simulation thread
	void Attach(base::SPHSimulation<2, Particles>* sim) override;

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
	/// Executes on the UI thread
	void BeginFullscreenDockspace();
	/// Executes on the UI thread
	void BuildInitialLayout();

private:
	Particles m_Particles;

	GLFWwindow* window = nullptr;
	std::atomic<bool> m_Running = false;
	std::thread m_ImguiThread;
	std::mutex m_Mutex;
	unsigned int m_DockspaceID = 0, m_DockIDLeft = 0, m_DockIDCenter = 0;

	bool m_Changed = false;
	bool m_RequestNewParticles = true;

	float m_SimFPS = 0, m_SimTrueFPS = 0;
	int m_SimFramesCounter = 0;
	stdc::nanoseconds m_SimTimeCounter = 0ns, m_SimTrueTimeCounter = 0ns;
	stdc::time_point<stdclock> m_SimFrameStart;
	float m_MaxVelocity = 2.5, m_MaxPressure = 10000;
	Command<2> m_Cmd;

	float m_SimTime = 0;
	base::SPHProfiling m_SimProfiling;
	base::SPHParams m_SimParams;
	std::string m_SimName;

	enum ColoringParam
	{
		PRESSURE, VELOCITY
	};
	ColoringParam m_ColoringParam = VELOCITY;
};

#include "ImGuiViewer_Impl.hpp"