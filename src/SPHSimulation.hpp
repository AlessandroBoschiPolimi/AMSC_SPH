#pragma once
#include <vector>
#include <cmath>
#include "Observer.hpp"
#include "Kernel.hpp"

#include "Command.hpp"
#include "Initializers/SimInitializer.hpp"
#include "Neighbors/NeighborFinder.hpp"




template <size_t D>
class SPHSimulation
{
public:
	static constexpr size_t size = D;
	using      idx_t = Particle<D>::idx_t;
	using      vec_t = Particle<D>::vec_t;

	struct Params
	{
		float RestDensity = 1000.0f;
		float Stiffness = 1e2f;
		float Viscosity = 1e-4f;
		float ViscosityRigid = 5e-3f;
		float TimeStep = 0.0006f;
		float SmoothingLength = 0.007f;
		float PressureTol = 1e-2f;
		float FinalTime = 10.0;
	};

	struct Profiling
	{
		stdc::nanoseconds Neighbors, Initialize, IterativePressure;
	};

	friend class NeighborFinder<D>;

public: // Simulation interface
	SPHSimulation(NeighborFinder<D>* nf);
	~SPHSimulation() {
		for (auto* o : m_Observers)
			o->Attach(nullptr);
	}

	void InitializeFluid(const SimInitializer<D>* init);
	void Start();



private: // Simulation functions
	/*void BuildYWall(float x, float begin, float end, float delta);
	void BuildXWall(float y, float begin, float end, float delta);*/

	void Step();


	// Functions handling the three parts of the scheme for all the particles
	void FindAllNeighbors();
	void Initialize();
	void IterativePressure();
	// Helper functions for the subsequent particles
	void ComputeBoundaryPsi(idx_t i);
	void ComputeDensity(idx_t i);
	void ComputePressure(idx_t i);
	void ComputeAccelerationPressure(idx_t i);
	void ComputeAccelerationViscosity(idx_t i);
	void UpdatePositionInitial(idx_t i);
	void UpdatePositionIteration(idx_t i);
	void UpdateVelocityInitial(idx_t i);
	void UpdateVelocityIteration(idx_t i);

	void EvaluateCommand(idx_t i);


private:
	std::vector<Observer<D>*> m_Observers;
	std::vector<Particle<D>> m_Particles;

	size_t m_Frame = 0;
	float m_Time = 0.0f;

	Profiling m_Profiling;

	Params m_Params;
	Command<D> m_Command;

	NeighborFinder<D>* m_NeighborFinder = nullptr;

	Kernel<D> W_Ker;


public: // Generic interface
	void AddObserver(Observer<D>* obs) {
		obs->Attach(this);
		m_Observers.push_back(obs);
	}

	void SetParams(const Params& params) { m_Params = params; }
	void SetRestDensity    (float val) { m_Params.RestDensity     = val; }
	void SetStiffness      (float val) { m_Params.Stiffness       = val; }
	void SetViscosity      (float val) { m_Params.Viscosity       = val; }
	void SetTimeStep       (float val) { m_Params.TimeStep        = val; }
	void SetSmoothingLength(float val) { m_Params.SmoothingLength = val; }
	Params GetParams() const { return m_Params; }
	float GetRestDensity    () const { return m_Params.RestDensity    ; }
	float GetStiffness      () const { return m_Params.Stiffness      ; }
	float GetViscosity      () const { return m_Params.Viscosity      ; }
	float GetTimeStep       () const { return m_Params.TimeStep       ; }
	float GetSmoothingLength() const { return m_Params.SmoothingLength; }

	const std::vector<Particle<D>>& GetParticles() const { return m_Particles; }

	float  GetTime()  const { return m_Time; }
	size_t GetFrame() const { return m_Frame; }
	Profiling GetProfiling() const { return m_Profiling; }

	void ApplyCommand(const Command<D>& cmd) { m_Command = cmd; }


private: // Generic functions
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


// TODO: m_Params.SmoothingLength isn't garbage only if m_Params is defined before W_Ker, abstract the default value
template<size_t D>
inline SPHSimulation<D>::SPHSimulation(NeighborFinder<D>* nf)
	: m_NeighborFinder(nf), W_Ker(m_Params.SmoothingLength/ 2.0f)
{ }


struct Time
{
	Time(stdc::nanoseconds* res) {
		Res = res;
		Start = stdclock::now();
	}
	~Time() {
		*Res = stdclock::now() - Start;
	}

	stdc::nanoseconds* Res = nullptr;
	stdc::time_point<stdclock> Start;
};

template <size_t D>
inline void SPHSimulation<D>::Start()
{
	while (m_Time < m_Params.FinalTime)
		Step();
}
template <size_t D>
inline void SPHSimulation<D>::Step()
{
	NotifyStartFrame();

	{
		Time time(&m_Profiling.Neighbors);
		m_NeighborFinder->InitializeFrame(this);
		FindAllNeighbors();
	}
	{
		Time time(&m_Profiling.Initialize);
		Initialize();
	}
	{
		Time time(&m_Profiling.IterativePressure);
		IterativePressure();
	}

	m_Time += m_Params.TimeStep;
	m_Frame++;
	NotifyEndFrame();
}



template <size_t D>
inline void SPHSimulation<D>::FindAllNeighbors()
{
	for (int i = 0; i < m_Particles.size(); i++)
		m_NeighborFinder->Find(i, m_Particles[i].Neighbors);
}

template <size_t D>
inline void SPHSimulation<D>::Initialize()
{
	/*
	 * Non-iterative part of the timestep
	 * In each timestep, set initial force due to Viscosity
	 * Apply this and gravity force to all the particles
	 * Additionally, we need to compute 'Mass' of boundary particles (ParticlePsi)
	 * In the first step, we also need to initialize density
	 */
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == SOLID)
			ComputeBoundaryPsi(i);
	}
	if (m_Frame == 0)
	{
		for (int i = 0; i < m_Particles.size(); i++)
		{
			if (m_Particles[i].Type == FLUID)
				ComputeDensity(i);
		}
	}
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == FLUID) {
			ComputeAccelerationViscosity(i);
		}
	}
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == FLUID) {
			UpdateVelocityInitial(i);
			UpdatePositionInitial(i);
		}
	}
}
template <size_t D>
inline void SPHSimulation<D>::IterativePressure()
{
	/*
	 * Use simple scheme with splitting
	 * After computing initial forces and moving particles, compute dansity and pressure
	   and move particles again.
	 */
	for (int i = 0; i < m_Particles.size(); i++) {
		if (m_Particles[i].Type == SOLID)
			continue;
		ComputeDensity(i);
		ComputePressure(i);
	}

	if (m_Command.Type != Command<D>::NONE)
		for (int i = 0; i < m_Particles.size(); i++)
			EvaluateCommand(i);

	for (int i = 0; i < m_Particles.size(); i++) {
		if (m_Particles[i].Type == SOLID)
			continue;
		ComputeAccelerationPressure(i);
		UpdateVelocityIteration(i);
		UpdatePositionIteration(i);
	}
}

template <size_t D>
inline void SPHSimulation<D>::ComputeBoundaryPsi(idx_t i)
{
	/*
	 * Computes 'mass' of the boundary particles used
	   to implement collisions
	 */
	float V = 0;
	for (auto& j : m_Particles[i].Neighbors)
	{
		if (m_Particles[j].Type == SOLID)
			V += W_Ker.GetValue(m_Particles[i], m_Particles[j]);
	}
	//Clamp the values in case the volume is too small
	m_Particles[i].BoundaryPsi = (V > 1.0f) ? m_Params.RestDensity / V : 0;
}
template <size_t D>
inline void SPHSimulation<D>::ComputeDensity(idx_t i)
{
	/*
	 * Simple function to compute density
	 * Takes initial value of density produced by itself
	 * For security, clamps the value in the end to avoid disappearing particles
	 */
	m_Particles[i].Density = m_Particles[i].Mass * W_Ker.GetValue(m_Particles[i], m_Particles[i]);
	for (auto& j : m_Particles[i].Neighbors) {
		float W_ij = W_Ker.GetValue(m_Particles[i], m_Particles[j]);
		if (m_Particles[j].Type == FLUID)
		{
			m_Particles[i].Density += m_Particles[j].Mass * W_ij;
		}
		//Boundary handling
		else
		{
			m_Particles[i].Density += m_Particles[j].BoundaryPsi * W_ij;
		}
	}
}

template <size_t D>
inline void SPHSimulation<D>::ComputePressure(idx_t i)
{
	/*
	 * (Andrew) Tait equation
	 * Stiffness constant is user defined
	 */
	m_Particles[i].Pressure = std::max(m_Params.Stiffness *
		(std::pow(m_Particles[i].Density /
			m_Params.RestDensity, 7.0f) - 1), 0.0f);
}

template <size_t D>
inline void SPHSimulation<D>::ComputeAccelerationViscosity(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [1]
	 * [2] has typos in that formula
	 */
	m_Particles[i].A_visc = vec_t{ 0, 0 };
	for (auto& j : m_Particles[i].Neighbors)
	{
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		vec_t v_ij = m_Particles[i].Velocity - m_Particles[j].Velocity;
		vec_t x_ij = m_Particles[i].Position - m_Particles[j].Position;
		float prod = Dot(x_ij, v_ij);
		float numerator = ((prod > 0) ? prod : 0);
		if (m_Particles[j].Type == FLUID)
		{
			numerator *= m_Particles[j].Mass;
			float denominator = (m_Particles[i].Density + m_Particles[j].Density) * (Dot(x_ij, x_ij) +
				0.01 * m_Params.SmoothingLength * m_Params.SmoothingLength);
			m_Particles[i].A_visc += (2 * m_Params.Viscosity * (numerator / denominator)) * DW_ij;
		}
		else
		{
			float denominator = 2 * m_Particles[i].Density * (Dot(x_ij, x_ij) +
				0.01 * m_Params.SmoothingLength * m_Params.SmoothingLength);
			m_Particles[i].A_visc += (m_Params.ViscosityRigid * (numerator / denominator)) *
				m_Particles[j].BoundaryPsi * DW_ij;
		}
	}
}

template <size_t D>
inline void SPHSimulation<D>::ComputeAccelerationPressure(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [2]
	 * Boundary handling according to [2]
	 */
	m_Particles[i].A_press = vec_t{ 0, 0 };
	for (auto& j : m_Particles[i].Neighbors)
	{
		float factor = (m_Particles[j].Type == SOLID) ? (m_Particles[j].BoundaryPsi / 2.0f) : m_Particles[j].Mass;
		idx_t z = (m_Particles[j].Type == SOLID) ? i : j;
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		m_Particles[i].A_press -=
			(m_Particles[i].Pressure / (m_Particles[i].Density * m_Particles[i].Density) +
				m_Particles[z].Pressure / (m_Particles[z].Density * m_Particles[z].Density)) *
			factor *
			DW_ij;
	}
}

template <size_t D>
inline void SPHSimulation<D>::UpdatePositionInitial(idx_t i)
{
	//Update position due to viscosity and gravity
	m_Particles[i].Position +=
		m_Params.TimeStep *
		m_Particles[i].Velocity;
}
template <size_t D>
inline void SPHSimulation<D>::UpdatePositionIteration(idx_t i)
{
	//Update position due to pressure
	m_Particles[i].Position +=
		m_Params.TimeStep *
		m_Params.TimeStep *
		m_Particles[i].A_press;
}
template <size_t D>
inline void SPHSimulation<D>::UpdateVelocityInitial(idx_t i)
{
	//Update velocity due to viscosity and gravity
	m_Particles[i].Velocity +=
		m_Params.TimeStep *
		(m_Particles[i].A_grav +
			m_Particles[i].A_visc);
}
template <size_t D>
inline void SPHSimulation<D>::UpdateVelocityIteration(idx_t i)
{
	//Update velocity due to pressure
	m_Particles[i].Velocity +=
		m_Params.TimeStep *
		m_Particles[i].A_press;
}


template <size_t D>
inline void SPHSimulation<D>::EvaluateCommand(idx_t i)
{
	float d = Norm(m_Particles[i].Position - m_Command.Position);
	if (d < m_Command.Radius)
	{
		float falloff = 1.0f - d / m_Command.Radius;
		m_Particles[i].Pressure += falloff * m_Command.Strength;
	}
}





template<size_t D>
inline void SPHSimulation<D>::InitializeFluid(const SimInitializer<D>* init)
{
	std::cout << "Initializing" << '\n';
	init->Init(m_Particles, m_Params.SmoothingLength);
	std::cout << "Particles: " << m_Particles.size() << '\n';
}
