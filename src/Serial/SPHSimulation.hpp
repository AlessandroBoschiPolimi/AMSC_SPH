#pragma once
#include "Neighbors/NeighborFinder.hpp"
#include "Base/SPHSimulation.hpp"


namespace serial
{

template <size_t D, ParticleSet<D> Particles>
class SPHSimulation : public base::SPHSimulation<D, Particles>
{
public:
	using idx_t = Particle<D>::idx_t;
	using vec_t = Particle<D>::vec_t;

	friend class NeighborFinder<D, Particles>;


	SPHSimulation(NeighborFinder<D, Particles>* nf);
	virtual ~SPHSimulation() override { }

	void Start() override;


	const Particles& GetParticles() const override { return m_Particles; }
	Particles& GetParticles() override { return m_Particles; }
	const std::vector<std::vector<idx_t>>& GetNeighbors() const override { return m_Neighbors; }
	std::vector<std::vector<idx_t>>& GetNeighbors() override { return m_Neighbors; }

	void InitializeFluid(const SimInitializer<D, Particles>* init) override;


protected:
	void Step() override;


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
	Particles m_Particles;
	std::vector<std::vector<idx_t>> m_Neighbors;

	NeighborFinder<D, Particles>* m_NeighborFinder = nullptr;

	Kernel<D> W_Ker;
};


/* References:
 * [1]https://cg.informatik.uni-freiburg.de/publications/2012_SIGGRAPH_rigidFluidCoupling.pdfa
 * [2]https://cg.informatik.uni-freiburg.de/course_notes/sim_10_sph.pdf
 */


template <size_t D, ParticleSet<D> Particles>
inline SPHSimulation<D, Particles>::SPHSimulation(NeighborFinder<D, Particles>* nf)
	: base::SPHSimulation<D, Particles>(), m_NeighborFinder(nf), W_Ker(this->m_Params.SmoothingLength / 2.0f)
{ }



template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Start()
{
	std::cout << "Start!\n";
	while (this->m_Time < this->m_Params.FinalTime)
		Step();
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Step()
{
	this->NotifyStartFrame();
	for (auto& obj : this->m_Objects)
		obj->OnFrameStart(this->m_Time);

	if (m_Neighbors.size() != m_Particles.Size())
	{
		m_Neighbors.clear();
		m_Neighbors.resize(m_Particles.Size());
	}

	{
		stdc::time_point<stdclock> start;
		start = stdclock::now();

		m_NeighborFinder->InitializeFrame(this);
		FindAllNeighbors();

		this->m_Profiling.Neighbors = stdclock::now() - start;
	}
	{
		stdc::time_point<stdclock> start;
		start = stdclock::now();

		Initialize();

		this->m_Profiling.Initialize = stdclock::now() - start;
	}
	{
		stdc::time_point<stdclock> start;
		start = stdclock::now();

		IterativePressure();

		this->m_Profiling.IterativePressure = stdclock::now() - start;
	}

	this->m_Time += this->m_Params.TimeStep;
	this->m_Frame++;
	this->NotifyEndFrame();
}



template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::FindAllNeighbors()
{
	for (int i = 0; i < m_Particles.Size(); i++)
	{
		if (m_Particles.Type(i) == SOLID && this->m_Frame > 0)
			continue;
		m_NeighborFinder->Find(i, m_Neighbors[i]);
	}
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Initialize()
{
	/*
	 * Non-iterative part of the timestep
	 * In each timestep, set initial force due to Viscosity
	 * Apply this and gravity force to all the particles
	 * Additionally, we need to compute 'Mass' of boundary particles (ParticlePsi)
	 * In the first step, we also need to initialize density
	 */
	if (this->m_Frame == 0)
	{
		for (int i = 0; i < m_Particles.Size(); i++)
		{
			if (m_Particles.Type(i) == FLUID)
				ComputeDensity(i);
			else
				ComputeBoundaryPsi(i);
		}
	}
	for (int i = 0; i < m_Particles.Size(); i++)
	{
		if (m_Particles.Type(i) == FLUID)
			ComputeAccelerationViscosity(i);
	}
	for (int i = 0; i < m_Particles.Size(); i++)
	{
		if (m_Particles.Type(i) == FLUID) {
			UpdateVelocityInitial(i);
			UpdatePositionInitial(i);
		}
	}
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::IterativePressure()
{
	/*
	 * Use simple scheme with splitting
	 * After computing initial forces and moving particles, compute dansity and pressure
	   and move particles again.
	 */

	for (int i = 0; i < m_Particles.Size(); i++)
	{
		if (m_Particles.Type(i) == SOLID)
			continue;
		ComputeDensity(i);
		ComputePressure(i);
	}
	if (this->m_Command.Type != Command<D>::NONE)
	{
		for (int i = 0; i < m_Particles.Size(); i++)
			EvaluateCommand(i);
	}
	for (int i = 0; i < m_Particles.Size(); i++)
	{
		if (m_Particles.Type(i) == SOLID)
			continue;
		ComputeAccelerationPressure(i);
	}
	for (int i = 0; i < m_Particles.Size(); i++)
	{
		if (m_Particles.Type(i) == SOLID)
			continue;
		UpdatePositionIteration(i);
		UpdateVelocityIteration(i);
	}
}


template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputeBoundaryPsi(idx_t i)
{
	/*
	 * Computes 'mass' of the boundary particles used
	 * to implement collisions
	 */
	float V = 0;
	for (auto& j : m_Neighbors[i])
	{
		if (m_Particles.Type(j) == SOLID)
			V += W_Ker.GetValue(m_Particles.Position(i), m_Particles.Position(j));
	}
	// Clamp the values in case the volume is too small
	m_Particles.SetBoundaryPsi(i, (V > 1.0f) ? this->m_Params.RestDensity / V : 0);
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputeDensity(idx_t i)
{
	/*
	 * Simple function to compute density
	 * Takes initial value of density produced by itself
	 * For security, clamps the value in the end to avoid disappearing particles
	 */
	
	auto pi = m_Particles.Position(i);

	float density = m_Particles.Mass(i) * W_Ker.GetValue(pi, pi);
	for (auto& j : m_Neighbors[i])
	{
		float W_ij = W_Ker.GetValue(pi, m_Particles.Position(j));
		if (m_Particles.Type(j) == FLUID)
			density += m_Particles.Mass(j) * W_ij;
		else // Boundary handling
			density += m_Particles.BoundaryPsi(j) * W_ij;
	}

	m_Particles.SetDensity(i, density);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputePressure(idx_t i)
{
	/*
	 * (Andrew) Tait equation
	 * Stiffness constant is user defined
	 */

	float pressure = std::max(
						this->m_Params.Stiffness * (std::pow(m_Particles.Density(i) / (this->m_Params.RestDensity), 7.0f) - 1),
						0.0f
					 );
	m_Particles.SetPressure(i, pressure);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputeAccelerationViscosity(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [1]
	 * [2] has typos in that formula
	 */

	auto pi = m_Particles.Position(i);
	auto vi = m_Particles.Velocity(i);
	auto di = m_Particles.Density(i);

	float sl2 = this->m_Params.SmoothingLength * this->m_Params.SmoothingLength;
	vec_t acc;
	if constexpr (D == 2)
		acc = vec_t{ 0, 0 };
	else
		acc = vec_t{ 0, 0, 0 };

	for (auto& j : m_Neighbors[i])
	{
		vec_t DW_ij = W_Ker.GetGradient(pi, m_Particles.Position(j));
		vec_t v_ij = vi - m_Particles.Velocity(j);
		vec_t x_ij = pi - m_Particles.Position(j);
		float prod = Dot(x_ij, v_ij);
		float numerator = ((prod > 0) ? prod : 0);

		if (m_Particles.Type(j) == FLUID)
		{
			numerator *= m_Particles.Mass(j);
			float denominator = (di + m_Particles.Density(j)) * (Dot(x_ij, x_ij) + 0.01 * sl2);
			acc += (2 * this->m_Params.Viscosity * (numerator / denominator)) * DW_ij;
		}
		else
		{
			float denominator = 2 * di * (Dot(x_ij, x_ij) + 0.01 * sl2);
			acc += (this->m_Params.ViscosityRigid * (numerator / denominator)) * m_Particles.BoundaryPsi(j) * DW_ij;
		}
	}

	m_Particles.SetAccelerationViscosity(i, acc);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputeAccelerationPressure(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [2]
	 * Boundary handling according to [2]
	 */

	auto pi = m_Particles.Position(i);
	auto pri = m_Particles.Pressure(i);
	auto di = m_Particles.Density(i);
	auto ddi = di * di;

	vec_t acc;
	if constexpr (D == 2)
		acc = vec_t{ 0, 0 };
	else
		acc = vec_t{ 0, 0, 0 };

	for (auto& j : m_Neighbors[i])
	{
		float factor = (m_Particles.Type(j) == SOLID) ? (m_Particles.BoundaryPsi(j) / 2.0f) : m_Particles.Mass(j);
		idx_t z = (m_Particles.Type(j) == SOLID) ? i : j;
		vec_t DW_ij = W_Ker.GetGradient(pi, m_Particles.Position(j));
		auto dz = m_Particles.Density(z);

		acc -= (pri / ddi + m_Particles.Pressure(z) / (dz * dz)) * factor * DW_ij;
	}

	m_Particles.SetAccelerationPressure(i, acc);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdatePositionInitial(idx_t i)
{
	// Update position due to viscosity and gravity
	auto pi = m_Particles.Position(i);
	m_Particles.SetPosition(i, pi + this->m_Params.TimeStep * m_Particles.Velocity(i));
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdatePositionIteration(idx_t i)
{
	// Update position due to pressure
	auto pi = m_Particles.Position(i);
	pi += this->m_Params.TimeStep * this->m_Params.TimeStep *
		  m_Particles.AccelerationPressure(i);
	m_Particles.SetPosition(i, pi);

	for (auto &obj: this->m_Objects)
		obj->Activate(i);
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdateVelocityInitial(idx_t i)
{
	// Update velocity due to viscosity and gravity
	auto vi = m_Particles.Velocity(i);
	vi += this->m_Params.TimeStep * (m_Particles.AccelerationGravity(i) + m_Particles.AccelerationViscosity(i));
	m_Particles.SetVelocity(i, vi);
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdateVelocityIteration(idx_t i)
{
	// Update velocity due to pressure
	auto vi = m_Particles.Velocity(i);
	vi += this->m_Params.TimeStep * m_Particles.AccelerationPressure(i);
	m_Particles.SetVelocity(i, vi);
}


template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::EvaluateCommand(idx_t i)
{
	float d = Norm(m_Particles.Position(i) - this->m_Command.Position);
	if (d < this->m_Command.Radius)
	{
		float falloff = 1.0f - d / this->m_Command.Radius;
		float pressure = m_Particles.Pressure(i) + falloff * this->m_Command.Strength;
		m_Particles.SetPressure(i, pressure);
	}
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::InitializeFluid(const SimInitializer<D, Particles>* init)
{
	std::cout << "Initializing " << this->m_Name << '\n';
	init->Init(m_Particles, this->m_Params.SmoothingLength, this->m_Objects);
	std::cout << "Particles: " << m_Particles.Size() << '\n';
}

}