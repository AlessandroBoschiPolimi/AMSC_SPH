#pragma once
#include <vector>
#include "Particle.hpp"

class SPHSimulation;

class Observer {
public:
	Observer() = default;
	virtual ~Observer() = default;

	/// Mainly used to read simulation parameters / timing
	/// Executes on the simulation thread
	virtual void OnStartFrame() {}
	/// Mainly used to alter simulation parameters / timing
	/// Executes on the simulation thread
	virtual void OnEndFrame() {}

	/// Executes on the simulation thread
	virtual void Attach(SPHSimulation* sim) { m_Sim = sim; }

protected:
	SPHSimulation* m_Sim = nullptr;
};
