#pragma once
#include "Neighbors/NeighborFinder.hpp"
#include "Base/SPHSimulation.hpp"



namespace openmp
{

template <size_t D>
class SPHSimulation : public base::SPHSimulation<D>
{
public:
	using idx_t = Particle<D>::idx_t;
	using vec_t = Particle<D>::vec_t;

	friend class NeighborFinder<D>;


	SPHSimulation(NeighborFinder<D>* nf);
	virtual ~SPHSimulation() override { }

	void Start() override;


	const std::vector<Particle<D>>& GetParticles() const override { return m_Particles; }
	std::vector<Particle<D>>& GetParticles() override { return m_Particles; }
	const std::vector<std::vector<idx_t>>& GetNeighbors() const override { return m_Neighbors; }
	std::vector<std::vector<idx_t>>& GetNeighbors() override { return m_Neighbors; }

	void InitializeFluid(const SimInitializer<D>* init) override;


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
	std::vector<Observer<D>*> m_Observers;
	std::vector<Particle<D>> m_Particles;
	std::vector<std::vector<idx_t>> m_Neighbors;

	NeighborFinder<D>* m_NeighborFinder = nullptr;

	Kernel<D> W_Ker;
};


/* References:
 * [1]https://cg.informatik.uni-freiburg.de/publications/2012_SIGGRAPH_rigidFluidCoupling.pdfa
 * [2]https://cg.informatik.uni-freiburg.de/course_notes/sim_10_sph.pdf
 */


template<size_t D>
inline SPHSimulation<D>::SPHSimulation(NeighborFinder<D>* nf) 
	: m_NeighborFinder(nf), W_Ker(this->m_Params.SmoothingLength / 2.0f)
{ }



template <size_t D>
inline void SPHSimulation<D>::Start()
{
	#pragma omp parallel
	{
		while (this->m_Time < this->m_Params.FinalTime)
			Step();
	}
}
template <size_t D>
inline void SPHSimulation<D>::Step()
{
	#pragma omp single
	{
		this->NotifyStartFrame();

		for (auto& obj : this->m_Objects)
			obj->Activate();

		if (m_Neighbors.size() != m_Particles.size())
		{
			m_Neighbors.clear();
			m_Neighbors.resize(m_Particles.size());
		}
	}

	{
		stdc::time_point<stdclock> start;
		#pragma omp master
		{ start = stdclock::now(); }
		#pragma omp barrier

		m_NeighborFinder->InitializeFrame(this);
		#pragma omp barrier
		FindAllNeighbors();

		#pragma omp barrier
		#pragma omp master
		{ this->m_Profiling.Neighbors = stdclock::now() - start; }
		#pragma omp barrier
	}
	{
		stdc::time_point<stdclock> start;
		#pragma omp master
		{ start = stdclock::now(); }
		#pragma omp barrier

		Initialize();

		#pragma omp master
		{ this->m_Profiling.Initialize = stdclock::now() - start; }
		#pragma omp barrier
	}
	{
		stdc::time_point<stdclock> start;
		#pragma omp master
		{ start = stdclock::now(); }
		#pragma omp barrier

		IterativePressure();

		#pragma omp master
		{ this->m_Profiling.IterativePressure = stdclock::now() - start; }
		#pragma omp barrier
	}

	#pragma omp single
	{
		this->m_Time += this->m_Params.TimeStep;
		this->m_Frame++;
		this->NotifyEndFrame();
	}
}



template <size_t D>
inline void SPHSimulation<D>::FindAllNeighbors()
{
	#pragma omp for
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == SOLID && this->m_Frame > 0)
			continue;
		this->m_NeighborFinder->Find(i, this->m_Neighbors[i]);
	}
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
	if (this->m_Frame == 0)
	{
		#pragma omp for
		for (int i = 0; i < m_Particles.size(); i++)
		{
			if (m_Particles[i].Type == FLUID)
				ComputeDensity(i);
			else
				ComputeBoundaryPsi(i);
		}
	}
	#pragma omp for
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == FLUID) {
			ComputeAccelerationViscosity(i);
		}
	}
	#pragma omp for
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

	#pragma omp for
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == SOLID)
			continue;
		ComputeDensity(i);
		ComputePressure(i);
	}
	if (this->m_Command.Type != Command<D>::NONE)
	{
		#pragma omp for
		for (int i = 0; i < m_Particles.size(); i++)
			EvaluateCommand(i);
	}
	#pragma omp for
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == SOLID)
			continue;
		ComputeAccelerationPressure(i);
	}
	#pragma omp for
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == SOLID)
			continue;
		UpdatePositionIteration(i);
		UpdateVelocityIteration(i);
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
	for (auto& j : m_Neighbors[i])
	{
		if (m_Particles[j].Type == SOLID)
			V += W_Ker.GetValue(m_Particles[i], m_Particles[j]);
	}
	//Clamp the values in case the volume is too small
	m_Particles[i].BoundaryPsi = (V > 1.0f) ? this->m_Params.RestDensity / V : 0;
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
	for (auto& j : m_Neighbors[i]) {
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
	m_Particles[i].Pressure = std::max(this->m_Params.Stiffness *
		(std::pow(m_Particles[i].Density /
			(this->m_Params.RestDensity), 7.0f) - 1), 0.0f);
}

template <size_t D>
inline void SPHSimulation<D>::ComputeAccelerationViscosity(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [1]
	 * [2] has typos in that formula
	 */
	m_Particles[i].A_visc = vec_t{ 0, 0 };
	for (auto& j : m_Neighbors[i])
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
				0.01 * this->m_Params.SmoothingLength * this->m_Params.SmoothingLength);
			m_Particles[i].A_visc += (2 * this->m_Params.Viscosity * (numerator / denominator)) * DW_ij;
		}
		else
		{
			float denominator = 2 * m_Particles[i].Density * (Dot(x_ij, x_ij) +
				0.01 * this->m_Params.SmoothingLength * this->m_Params.SmoothingLength);
			m_Particles[i].A_visc += (this->m_Params.ViscosityRigid * (numerator / denominator)) *
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
	for (auto& j : m_Neighbors[i])
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
		this->m_Params.TimeStep *
		m_Particles[i].Velocity;
}
template <size_t D>
inline void SPHSimulation<D>::UpdatePositionIteration(idx_t i)
{
	//Update position due to pressure
	m_Particles[i].Position +=
		this->m_Params.TimeStep *
		this->m_Params.TimeStep *
		m_Particles[i].A_press;

	for (auto& obj : this->m_Objects)
		obj->Activate(i);
}
template <size_t D>
inline void SPHSimulation<D>::UpdateVelocityInitial(idx_t i)
{
	//Update velocity due to viscosity and gravity
	m_Particles[i].Velocity +=
		this->m_Params.TimeStep *
		(m_Particles[i].A_grav +
			m_Particles[i].A_visc);
}
template <size_t D>
inline void SPHSimulation<D>::UpdateVelocityIteration(idx_t i)
{
	//Update velocity due to pressure
	m_Particles[i].Velocity +=
		this->m_Params.TimeStep *
		m_Particles[i].A_press;
}


template <size_t D>
inline void SPHSimulation<D>::EvaluateCommand(idx_t i)
{
	float d = Norm(m_Particles[i].Position - this->m_Command.Position);
	if (d < this->m_Command.Radius)
	{
		float falloff = 1.0f - d / this->m_Command.Radius;
		m_Particles[i].Pressure += falloff * this->m_Command.Strength;
	}
}


template<size_t D>
inline void SPHSimulation<D>::InitializeFluid(const SimInitializer<D>* init)
{
	std::cout << "Initializing" << '\n';
	init->Init(m_Particles, this->m_Params.SmoothingLength, this->m_Objects);
	std::cout << "Particles: " << m_Particles.size() << '\n';
}


}