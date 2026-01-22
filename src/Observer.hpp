#pragma once
#include <vector>
#include "Particle.hpp"

class SPHSimulation;

class Observer {
public:
	Observer() = default;
	virtual ~Observer() = default;

	/// Mainly used to read simulation parameters / timing
	virtual void OnStartFrame() {}
	/// Mainly used to alter simulation parameters / timing
	virtual void OnEndFrame() {}

	virtual void Attach(SPHSimulation* sim) { m_Sim = sim; }

protected:
	SPHSimulation* m_Sim = nullptr;
};
