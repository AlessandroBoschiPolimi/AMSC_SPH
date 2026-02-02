#pragma once
#include <vector>
#include "Particle.hpp"

template <size_t D>
class SPHSimulation;

template <size_t D>
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
	virtual void Attach(SPHSimulation<D>* sim) { m_Sim = sim; }

protected:
	SPHSimulation<D>* m_Sim = nullptr;
};
