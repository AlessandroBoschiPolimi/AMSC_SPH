#pragma once
#include <vector>
#include "Particle.hpp"

namespace base
{
	template <size_t D, ParticleSet<D> Particles>
	class SPHSimulation;
}

template <size_t D, ParticleSet<D> Particles>
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
	virtual void Attach(base::SPHSimulation<D, Particles>* sim) { m_Sim = sim; }

protected:
	base::SPHSimulation<D, Particles>* m_Sim = nullptr;
};
