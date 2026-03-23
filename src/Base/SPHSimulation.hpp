#pragma once
#include <memory>
#include <vector>
#include <cmath>
#include <omp.h>

#include "Observer.hpp"
#include "Command.hpp"
#include "Kernel.hpp"

#include "Initializers/SimInitializer.hpp"
#include "Neighbors/NeighborFinder.hpp"

#include "Objects/Sink.hpp"
#include "Objects/Source.hpp"


namespace base
{


struct SPHParams
{
	float RestDensity = 1000.0f;
	float Stiffness = 3e2f;
	float Viscosity = 1e-5f;
	float ViscosityRigid = 5e-3f;
	float TimeStep = 0.0002f;
	float SmoothingLength = 0.007f;
	float PressureTol = 1e-2f;
	float FinalTime = 100.0;
};

struct SPHProfiling
{
	stdc::nanoseconds Neighbors = 0ns, Initialize = 0ns, IterativePressure = 0ns;
};

template <size_t D, ParticleSet<D> Particles>
class SPHSimulation
{
public:
	static constexpr size_t size = D;
	using idx_t = Particle<D>::idx_t;
	

public: // Simulation interface
	SPHSimulation() = default;
	virtual ~SPHSimulation() {
		for (auto* o : m_Observers)
			o->Attach(nullptr);
	}

	virtual void InitializeFluid(const SimInitializer<D, Particles>* init) = 0;
	virtual void Start() = 0;

protected: // Simulation functions

	virtual void Step() = 0;


protected:
	SPHParams m_Params;

	SPHProfiling m_Profiling;
	u64 m_Frame = 0;
	float m_Time = 0.0f;
	
	Command<D> m_Command;
	std::vector<std::unique_ptr<Object<D, Particles>>> m_Objects;

	std::vector<Observer<D, Particles>*> m_Observers;

	std::string m_Name = "Simulation";

public: // Generic interface
	void AddObserver(Observer<D, Particles>* obs) {
		obs->Attach(this);
		m_Observers.push_back(obs);
	}

	void SetName(const std::string& name) { m_Name = name; }
	std::string GetName() const { return m_Name; }

	void SetParams(const SPHParams& params) { m_Params = params; }
	void SetRestDensity     (float val)     { m_Params.RestDensity     = val; }
	void SetStiffness       (float val)     { m_Params.Stiffness       = val; }
	void SetViscosity       (float val)     { m_Params.Viscosity       = val; }
	void SetTimeStep        (float val)     { m_Params.TimeStep        = val; }
	void SetSmoothingLength (float val)     { m_Params.SmoothingLength = val; }
	void SetFinalTime		(float val)		{ m_Params.FinalTime       = val; }
	SPHParams GetParams     () const        { return m_Params; }
	float GetRestDensity    () const        { return m_Params.RestDensity    ; }
	float GetStiffness      () const        { return m_Params.Stiffness      ; }
	float GetViscosity      () const        { return m_Params.Viscosity      ; }
	float GetTimeStep       () const        { return m_Params.TimeStep       ; }
	float GetSmoothingLength() const        { return m_Params.SmoothingLength; }
	float GetFinalTime		() const { return m_Params.FinalTime; }

	virtual const Particles& GetParticles() const = 0;
	virtual Particles& GetParticles() = 0;
	virtual const std::vector<std::vector<idx_t>>& GetNeighbors() const = 0;
	virtual std::vector<std::vector<idx_t>>& GetNeighbors() = 0;

	float GetTime()  const { return m_Time;  }
	u64   GetFrame() const { return m_Frame; }
	SPHProfiling GetProfiling() const { return m_Profiling; }

	void ApplyCommand(const Command<D>& cmd) { m_Command = cmd; }


protected: // Generic functions
	void NotifyStartFrame() {
		for (auto* o : m_Observers)
			o->OnStartFrame();
	}
	void NotifyEndFrame() {
		for (auto* o : m_Observers)
			o->OnEndFrame();
	}
};


/* References:
 * [1]https://cg.informatik.uni-freiburg.de/publications/2012_SIGGRAPH_rigidFluidCoupling.pdfa
 * [2]https://cg.informatik.uni-freiburg.de/course_notes/sim_10_sph.pdf
 */
}
