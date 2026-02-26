#pragma once
#include "Neighbors/NeighborFinder.hpp"
#include "Base/SPHSimulation.hpp"


namespace mpi
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
	Particles& GetParticlesLocal() { return m_Particles_local; }
	Particles& GetParticlesGhost() { return m_Particles_ghost; }
	const std::vector<std::vector<idx_t>>& GetNeighbors() const override { return m_Neighbors; }
	std::vector<std::vector<idx_t>>& GetNeighbors() override { return m_Neighbors; }

	void InitializeFluid(const SimInitializer<D, Particles>* init) override;
	void SetRank(const int mpi_rank_, const MPI_Comm mpi_comm_, const int mpi_size_){mpi_rank = mpi_rank_; mpi_comm = mpi_comm_; mpi_size = mpi_size_;}

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
	//MPI related functions
	void FindBounds();
	void SplitParticles();
	void GatherParticles();
	void ExchangeParticlesLeft();
	void ExchangeParticlesRight();
	void ExchangeSizes();
	void PickLocalGhost(idx_t i, idx_t &ind, Particles* &particle);

private:
	Particles m_Particles;
	std::vector<std::vector<idx_t>> m_Neighbors;
	NeighborFinder<D, Particles>* m_NeighborFinder = nullptr;

	Kernel<D> W_Ker;
	//MPI variables
	int mpi_rank;
	MPI_Comm mpi_comm;
	int mpi_size;
	float max_x, min_x;
	float left_bound, right_bound;
	Particles m_Particles_local;
	Particles m_Particles_ghost;
	std::vector<int> size_local;
};


/* References:
 * [1]https://cg.informatik.uni-freiburg.de/publications/2012_SIGGRAPH_rigidFluidCoupling.pdfa
 * [2]https://cg.informatik.uni-freiburg.de/course_notes/sim_10_sph.pdf
 */


 // TODO: m_Params.SmoothingLength isn't garbage only if m_Params is defined before W_Ker, abstract the default value

template <size_t D, ParticleSet<D> Particles>
inline SPHSimulation<D, Particles>::SPHSimulation(NeighborFinder<D, Particles>* nf)
	: base::SPHSimulation<D, Particles>(), m_NeighborFinder(nf), W_Ker(this->m_Params.SmoothingLength / 2.0f)
{ }



template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Start()
{
	FindBounds();
	SplitParticles();
	while (this->m_Time < this->m_Params.FinalTime)
		Step();
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Step()
{
	if (mpi_rank == 0)
		this->NotifyStartFrame();
	m_Neighbors.clear();
	m_Neighbors.resize(m_Particles_local.Size());
	
	{
		stdc::time_point<stdclock> start;
		{ start = stdclock::now(); }
	
		m_NeighborFinder->InitializeFrame(this);
		FindAllNeighbors();
	
		{ this->m_Profiling.Neighbors = stdclock::now() - start; }
	}
	{
		stdc::time_point<stdclock> start;
		{ start = stdclock::now(); }
	
		Initialize();
	
		{ this->m_Profiling.Initialize = stdclock::now() - start; }
	}
	{
		stdc::time_point<stdclock> start;
		{ start = stdclock::now(); }
	
		IterativePressure();
	
		{ this->m_Profiling.IterativePressure = stdclock::now() - start; }
	}
	this->m_Time += this->m_Params.TimeStep;
	MPI_Barrier(mpi_comm);
	GatherParticles();
	m_Particles_ghost.Clear();
	ExchangeParticlesLeft();
	ExchangeParticlesRight();
	ExchangeSizes();
	this->m_Frame++;
	if (mpi_rank == 0)
		this->NotifyEndFrame();
}



template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::FindAllNeighbors()
{
	for (int i = 0; i < m_Particles_local.Size(); i++)
	{
		if (m_Particles_local.Type(i) == SOLID && this->m_Frame > 0)
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
		for (int i = 0; i < m_Particles_local.Size(); i++)
		{
			if (m_Particles_local.Type(i) == FLUID)
				ComputeDensity(i);
			else
				ComputeBoundaryPsi(i);
		}
	}
	for (int i = 0; i < m_Particles_local.Size(); i++)
	{
		if (m_Particles_local.Type(i) == FLUID)
			ComputeAccelerationViscosity(i);
	}
	for (int i = 0; i < m_Particles_local.Size(); i++)
	{
		if (m_Particles_local.Type(i) == FLUID) {
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

	for (int i = 0; i < m_Particles_local.Size(); i++)
	{
		if (m_Particles_local.Type(i) == SOLID)
			continue;
		ComputeDensity(i);
		ComputePressure(i);
	}
	if (this->m_Command.Type != Command<D>::NONE)
	{
		for (int i = 0; i < m_Particles_local.Size(); i++)
			EvaluateCommand(i);
	}
	for (int i = 0; i < m_Particles_local.Size(); i++)
	{
		if (m_Particles_local.Type(i) == SOLID)
			continue;
		ComputeAccelerationPressure(i);
	}
	for (int i = 0; i < m_Particles_local.Size(); i++)
	{
		if (m_Particles_local.Type(i) == SOLID)
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
		Particles* particle;
		idx_t ind;
		PickLocalGhost(j, ind, particle);
		if (particle->Type(ind) == SOLID)
			V += W_Ker.GetValue(m_Particles_local.Position(i), particle->Position(ind));
	}
	// Clamp the values in case the volume is too small
	m_Particles_local.SetBoundaryPsi(i, (V > 1.0f) ? this->m_Params.RestDensity / V : 0);
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputeDensity(idx_t i)
{
	/*
	 * Simple function to compute density
	 * Takes initial value of density produced by itself
	 * For security, clamps the value in the end to avoid disappearing particles
	 */
	
	auto pi = m_Particles_local.Position(i);

	float density = m_Particles_local.Mass(i) * W_Ker.GetValue(pi, pi);
	for (auto& j : m_Neighbors[i])
	{
		Particles* particle;
		idx_t ind;
		PickLocalGhost(j, ind, particle);
		float W_ij = W_Ker.GetValue(pi, particle->Position(ind));
		if (particle->Type(ind) == FLUID)
			density += particle->Mass(ind) * W_ij;
		else // Boundary handling
			density += particle->BoundaryPsi(ind) * W_ij;
	}

	m_Particles_local.SetDensity(i, density);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputePressure(idx_t i)
{
	/*
	 * (Andrew) Tait equation
	 * Stiffness constant is user defined
	 */

	float pressure = std::max(
						this->m_Params.Stiffness * (std::pow(m_Particles_local.Density(i) / (this->m_Params.RestDensity), 7.0f) - 1),
						0.0f
					 );
	m_Particles_local.SetPressure(i, pressure);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputeAccelerationViscosity(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [1]
	 * [2] has typos in that formula
	 */

	auto pi = m_Particles_local.Position(i);
	auto vi = m_Particles_local.Velocity(i);
	auto di = m_Particles_local.Density(i);

	float sl2 = this->m_Params.SmoothingLength * this->m_Params.SmoothingLength;
	vec_t acc = vec_t{ 0, 0 };

	for (auto& j : m_Neighbors[i])
	{
		Particles* particle;
		idx_t ind;
		PickLocalGhost(j, ind, particle);
		vec_t DW_ij = W_Ker.GetGradient(pi, particle->Position(ind));
		vec_t v_ij = vi - particle->Velocity(ind);
		vec_t x_ij = pi - particle->Position(ind);
		float prod = Dot(x_ij, v_ij);
		float numerator = ((prod > 0) ? prod : 0);

		if (particle->Type(ind) == FLUID)
		{
			numerator *= particle->Mass(ind);
			float denominator = (di + particle->Density(ind)) * (Dot(x_ij, x_ij) + 0.01 * sl2);
			acc += (2 * this->m_Params.Viscosity * (numerator / denominator)) * DW_ij;
		}
		else
		{
			float denominator = 2 * di * (Dot(x_ij, x_ij) + 0.01 * sl2);
			acc += (this->m_Params.ViscosityRigid * (numerator / denominator)) * particle->BoundaryPsi(ind) * DW_ij;
		}
	}

	m_Particles_local.SetAccelerationViscosity(i, acc);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ComputeAccelerationPressure(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [2]
	 * Boundary handling according to [2]
	 */

	auto pi = m_Particles_local.Position(i);
	auto pri = m_Particles_local.Pressure(i);
	auto di = m_Particles_local.Density(i);
	auto ddi = di * di;

	vec_t acc = vec_t{ 0, 0 };

	for (auto& j : m_Neighbors[i])
	{
		Particles* particle;
		idx_t ind;
		PickLocalGhost(j, ind, particle);
		float factor = (particle->Type(ind) == SOLID) ? (particle->BoundaryPsi(ind) / 2.0f) : particle->Mass(ind);
		vec_t DW_ij = W_Ker.GetGradient(pi, particle->Position(ind));
		idx_t z = (particle->Type(ind) == SOLID) ? i : ind;
		particle = (particle->Type(ind) == SOLID) ? &m_Particles_local : particle;
		auto dz = particle->Density(z);

		acc -= (pri / ddi + particle->Pressure(z) / (dz * dz)) * factor * DW_ij;
	}

	m_Particles_local.SetAccelerationPressure(i, acc);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdatePositionInitial(idx_t i)
{
	// Update position due to viscosity and gravity
	auto pi = m_Particles_local.Position(i);
	m_Particles_local.SetPosition(i, pi + this->m_Params.TimeStep * m_Particles_local.Velocity(i));
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdatePositionIteration(idx_t i)
{
	// Update position due to pressure
	auto pi = m_Particles_local.Position(i);
	pi += this->m_Params.TimeStep * this->m_Params.TimeStep *
		  m_Particles_local.AccelerationPressure(i);
	m_Particles_local.SetPosition(i, pi);
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdateVelocityInitial(idx_t i)
{
	// Update velocity due to viscosity and gravity
	auto vi = m_Particles_local.Velocity(i);
	vi += this->m_Params.TimeStep * (m_Particles_local.AccelerationGravity(i) + m_Particles_local.AccelerationViscosity(i));
	m_Particles_local.SetVelocity(i, vi);
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::UpdateVelocityIteration(idx_t i)
{
	// Update velocity due to pressure
	auto vi = m_Particles_local.Velocity(i);
	vi += this->m_Params.TimeStep * m_Particles_local.AccelerationPressure(i);
	m_Particles_local.SetVelocity(i, vi);
}


template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::EvaluateCommand(idx_t i)
{
	float d = Norm(m_Particles_local.Position(i) - this->m_Command.Position);
	if (d < this->m_Command.Radius)
	{
		float falloff = 1.0f - d / this->m_Command.Radius;
		float pressure = m_Particles_local.Pressure(i) + falloff * this->m_Command.Strength;
		m_Particles_local.SetPressure(i, pressure);
	}
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::InitializeFluid(const SimInitializer<D, Particles>* init)
{
	std::cout << "Initializing" << '\n';
	init->Init(m_Particles, this->m_Params.SmoothingLength, this->m_Objects);
	std::cout << "Particles: " << m_Particles.Size() << '\n';
}

/*
************************
*MPI FUNCTIONS
************************
*/
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::FindBounds()
{
	if (mpi_rank == 0)
	{
		min_x = m_Particles.PositionX(0);
		max_x = min_x;
		for (idx_t i = 1; i < m_Particles.Size(); ++i)
		{
			float pos_x = m_Particles.PositionX(i);
			if ( pos_x > max_x)
				max_x = pos_x;
			if (pos_x < min_x)
				min_x = pos_x;
		}
		float chunk = (max_x - min_x) / static_cast<float>(mpi_size);
		if (chunk < 2 * this->m_Params.SmoothingLength)
		{
			std::cerr << "Smoothing Length is greater than size of a section, unable to continue!\n";
        		MPI_Abort(MPI_COMM_WORLD, 1);
		}
		std::vector<float> bounds;
		bounds.reserve(mpi_size);
		for (int i = 0; i < mpi_size; i++)
			bounds.emplace_back(min_x + chunk * i);
		bounds.emplace_back(max_x);
		for (int i = 0; i < mpi_size; i++)
		{
		MPI_Send(&bounds[i], 1 , MPI_FLOAT, i,
		0,  mpi_comm);

		MPI_Send(&bounds[i+1], 1 , MPI_FLOAT, i,
		1,  mpi_comm);
		}
	}
	MPI_Recv(&left_bound , 1 , MPI_FLOAT , 0 ,
	0, mpi_comm, MPI_STATUS_IGNORE);
	
	MPI_Recv(&right_bound , 1 , MPI_FLOAT , 0 ,
	1, mpi_comm, MPI_STATUS_IGNORE);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::SplitParticles()
{
	int size, size_ghost_left, size_ghost_right;
	size_local.resize(mpi_size);
	if (mpi_rank == 0)
	{
		float chunk = (max_x - min_x) / static_cast<float>(mpi_size);
		std::vector<std::vector<Particle<D>>> division;
		std::vector<std::vector<Particle<D>>> division_ghost_right;
		std::vector<std::vector<Particle<D>>> division_ghost_left;
		division.resize(mpi_size);
		division_ghost_left.resize(mpi_size);
		division_ghost_right.resize(mpi_size);
		for (idx_t i = 0; i < m_Particles.Size(); ++i)
		{
			float x = m_Particles.PositionX(i);
			int sector = std::clamp(static_cast<int>(std::floor((x - min_x ) / chunk)),
				       	0, mpi_size - 1);
			int sector_left = std::max(0, 
				static_cast<int>(std::floor( (x - this->m_Params.SmoothingLength - min_x) / chunk)));
			int sector_right = std::min(mpi_size - 1, 
			   static_cast<int>(std::floor( (x + this->m_Params.SmoothingLength - min_x) / chunk)));
			division[sector].push_back(m_Particles.GetParticle(i));
			if (sector_left != sector)
				division_ghost_left[sector_left].push_back(m_Particles.GetParticle(i));
			if (sector_right != sector)
				division_ghost_right[sector_right].push_back(m_Particles.GetParticle(i));
		}
		for (int i = 1; i < mpi_size; i++)
		{
			size = division[i].size();
			size_local[i] = size * sizeof(Particle<D>);
			size_ghost_right = division_ghost_right[i].size();
			size_ghost_left = division_ghost_left[i].size();
			//Send local particals to each process except for rank 0 -> handle locally
			MPI_Send(&size, 1 , MPI_INT, i,
			0,  mpi_comm);
			MPI_Send(division[i].data(), size * sizeof(Particle<D>) , MPI_BYTE, i,
			1,  mpi_comm);

			//Block to send the ghost elements
//			{
//				// (1) The left ghosts of 1 handled locally, instead send right ghosts from 0 to 1
//				if (i == 1)
//				{
//					int size_tmp = division_ghost_right[0].size();
//					MPI_Send(&size_tmp, 1 , MPI_INT, 1,
//					4,  mpi_comm);
//					MPI_Send(division_ghost_right[0].data(), division_ghost_right[0].size() * sizeof(Particle<D>) , MPI_BYTE, 1,
//					5,  mpi_comm);
//				}
//				//(2) Send left ghosts of each process
//				else
//				{
//					MPI_Send(&size_ghost_left, 1 , MPI_INT, i - 1,
//					2,  mpi_comm);
//					MPI_Send(division_ghost_left[i].data(), size_ghost_left * sizeof(Particle<D>) , MPI_BYTE, i - 1,
//					3,  mpi_comm);
//				}
//				//(3) Send left ghosts of all processes except for left most one
//				if (i != mpi_size - 1)
//				{
//					MPI_Send(&size_ghost_right, 1 , MPI_INT, i + 1,
//					4,  mpi_comm);
//					MPI_Send(division_ghost_right[i].data(), size_ghost_right * sizeof(Particle<D>) , MPI_BYTE, i + 1,
//					5,  mpi_comm);
//				}
//			}
		}
//		//Handle process 0 locally
		size_local[0] = division[0].size() * sizeof(Particle<D>);
		m_Particles_local.ParticlesVector = division[0];
//		m_Particles_ghost.ParticlesVector = division_ghost_left[1];
	}
//	//Block of Receiving messages of particles
	if (mpi_rank != 0)
	{
		//Each process receives their local particles
		{
			MPI_Recv(&size , 1 , MPI_INT , 0 ,
			0, mpi_comm, MPI_STATUS_IGNORE);
		
			m_Particles_local.ParticlesVector.resize(size);
		
			MPI_Recv(m_Particles_local.ParticlesVector.data(), size * sizeof(Particle<D>) , MPI_BYTE , 0 ,
			1, mpi_comm, MPI_STATUS_IGNORE);
		}
//
//		//Each process receives their ghost particles
//		{
//			//Receive block of (1) and (2)
//			{
//				MPI_Recv(&size_ghost_left , 1 , MPI_INT , 0 ,
//				2, mpi_comm, MPI_STATUS_IGNORE);
//
//				m_Particles_ghost.ParticlesVector.resize(size_ghost_left);
//
//				MPI_Recv(m_Particles_ghost.ParticlesVector.data(), size_ghost_left * sizeof(Particle<D>) , MPI_BYTE , 0,
//				3, mpi_comm, MPI_STATUS_IGNORE);
//			}
//			//Receive block of (3)
//			if (mpi_rank != mpi_size - 1)
//			{
//				MPI_Recv(&size_ghost_right , 1 , MPI_INT , 0 ,
//				4, mpi_comm, MPI_STATUS_IGNORE);
//
//				int old_size = m_Particles_ghost.Size();
//				m_Particles_ghost.ParticlesVector.resize(size_ghost_right + old_size);
//
//				MPI_Recv(m_Particles_ghost.ParticlesVector.data() + old_size, size_ghost_right * sizeof(Particle<D>) , MPI_BYTE , 0 ,
//				5, mpi_comm, MPI_STATUS_IGNORE);
//			}
//		}
	}
	MPI_Bcast(size_local.data() , mpi_size , MPI_INT ,
	0 ,mpi_comm) ;
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::GatherParticles()
{
	int size, size_ghost;
	if (mpi_rank == 0)
	{
		size = m_Particles.Size();
		m_Particles.ParticlesVector.clear();
		m_Particles.ParticlesVector.resize(size);
	}
	std::vector<int> size_sum;
	size_sum.resize(mpi_size);
	size_sum[0] = 0;
	for (int i = 1; i < mpi_size; ++i)
		size_sum[i] = size_sum[i - 1] + size_local[i - 1];
	MPI_Gatherv(m_Particles_local.ParticlesVector.data(),
	m_Particles_local.Size() * sizeof(Particle<D>),
	MPI_BYTE,
       	m_Particles.ParticlesVector.data(),
       	size_local.data(),
	size_sum.data(),
       	MPI_BYTE ,
	0,
       	mpi_comm);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ExchangeParticlesLeft()
{
	//Exchange messages with the left neighbour
	if (mpi_rank != 0)
	{
		std::vector<idx_t> to_erase;
		std::vector<Particle<D>> left_ghost, left_owned;
		for (idx_t i = 0; i < m_Particles_local.Size(); ++i)
		{
			float pos_x = m_Particles_local.PositionX(i);
			float h = this->m_Params.SmoothingLength;
			if (pos_x < left_bound + h && pos_x > left_bound)
				left_ghost.push_back(m_Particles_local.GetParticle(i));
			if (pos_x <= left_bound && pos_x > left_bound - h)
			{
				left_owned.push_back(m_Particles_local.GetParticle(i));
				//m_Particles_ghost.PushBack(m_Particles_local.GetParticle(i));
				to_erase.push_back(i);
			}
			if (pos_x < left_bound - h)
			{
				left_owned.push_back(m_Particles_local.GetParticle(i));
				to_erase.push_back(i);
			}
		}
		std::sort(to_erase.rbegin(), to_erase.rend());
		for (int idx : to_erase)
		{
    			m_Particles_local.ParticlesVector.erase(m_Particles_local.ParticlesVector.begin() + idx);
		}
		int size_owned = left_owned.size();
		int size_ghost = left_ghost.size();
		int receiver = mpi_rank - 1;
		MPI_Send(&size_owned, 1 , MPI_INT, receiver,
		0,  mpi_comm);
		MPI_Send(&size_ghost, 1 , MPI_INT, receiver,
		1,  mpi_comm);
		MPI_Send(left_owned.data(), size_owned * sizeof(Particle<D>) , MPI_BYTE, receiver,
		2,  mpi_comm);
		MPI_Send(left_ghost.data(), size_ghost * sizeof(Particle<D>) , MPI_BYTE, receiver,
		3,  mpi_comm);
	}
	if (mpi_rank != mpi_size - 1)
	{
		int size_owned, size_ghost;
		int sender = mpi_rank + 1;
		
		MPI_Recv(&size_owned, 1 , MPI_INT , sender ,
		0, mpi_comm, MPI_STATUS_IGNORE);
		
		int old_size_local = m_Particles_local.Size();
		m_Particles_local.ParticlesVector.resize(old_size_local + size_owned);

		MPI_Recv(&size_ghost, 1 , MPI_INT , sender ,
		1, mpi_comm, MPI_STATUS_IGNORE);
		
		int old_size_ghost = m_Particles_ghost.Size();
		m_Particles_ghost.ParticlesVector.resize(old_size_ghost + size_ghost);

		MPI_Recv(m_Particles_local.ParticlesVector.data() + old_size_local,
		size_owned * sizeof(Particle<D>) , MPI_BYTE , sender ,
		2, mpi_comm, MPI_STATUS_IGNORE);
		
		MPI_Recv(m_Particles_ghost.ParticlesVector.data() + old_size_ghost,
		size_ghost * sizeof(Particle<D>) , MPI_BYTE , sender ,
		3, mpi_comm, MPI_STATUS_IGNORE);

	}

}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ExchangeParticlesRight()
{
	//Exchange messages with the left neighbour
	if (mpi_rank != mpi_size - 1)
	{
		std::vector<idx_t> to_erase;
		std::vector<Particle<D>> right_ghost, right_owned;
		for (idx_t i = 0; i < m_Particles_local.Size(); ++i)
		{
			float pos_x = m_Particles_local.PositionX(i);
			float h = this->m_Params.SmoothingLength;
			if (pos_x > right_bound - h && pos_x <= right_bound)
				right_ghost.push_back(m_Particles_local.GetParticle(i));
			if (pos_x > right_bound && pos_x < right_bound + h)
			{
				right_owned.push_back(m_Particles_local.GetParticle(i));
				//m_Particles_ghost.PushBack(m_Particles_local.GetParticle(i));
				to_erase.push_back(i);
			}
			if (pos_x > right_bound + h)
			{
				right_owned.push_back(m_Particles_local.GetParticle(i));
				to_erase.push_back(i);
			}
		}
		std::sort(to_erase.rbegin(), to_erase.rend());
		for (int idx : to_erase)
		{
    			m_Particles_local.ParticlesVector.erase(m_Particles_local.ParticlesVector.begin() + idx);
		}
		int size_owned = right_owned.size();
		int size_ghost = right_ghost.size();
		int receiver = mpi_rank + 1;
		MPI_Send(&size_owned, 1 , MPI_INT, receiver,
		0,  mpi_comm);
		MPI_Send(&size_ghost, 1 , MPI_INT, receiver,
		1,  mpi_comm);
		MPI_Send(right_owned.data(), size_owned * sizeof(Particle<D>) , MPI_BYTE, receiver,
		2,  mpi_comm);
		MPI_Send(right_ghost.data(), size_ghost * sizeof(Particle<D>) , MPI_BYTE, receiver,
		3,  mpi_comm);
	}
	if (mpi_rank != 0)
	{
		int size_owned, size_ghost;
		int sender = mpi_rank - 1;
		
		MPI_Recv(&size_owned, 1 , MPI_INT , sender ,
		0, mpi_comm, MPI_STATUS_IGNORE);
		
		int old_size_local = m_Particles_local.Size();
		m_Particles_local.ParticlesVector.resize(old_size_local + size_owned);

		MPI_Recv(&size_ghost, 1 , MPI_INT , sender ,
		1, mpi_comm, MPI_STATUS_IGNORE);
		
		int old_size_ghost = m_Particles_ghost.Size();
		m_Particles_ghost.ParticlesVector.resize(old_size_ghost + size_ghost);

		MPI_Recv(m_Particles_local.ParticlesVector.data() + old_size_local,
		size_owned * sizeof(Particle<D>) , MPI_BYTE , sender ,
		2, mpi_comm, MPI_STATUS_IGNORE);
		
		
		MPI_Recv(m_Particles_ghost.ParticlesVector.data() + old_size_ghost,
		size_ghost * sizeof(Particle<D>) , MPI_BYTE , sender ,
		3, mpi_comm, MPI_STATUS_IGNORE);
	}

}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::ExchangeSizes()
{
	int my_size = m_Particles_local.Size() * sizeof(Particle<D>);
	MPI_Allgather(
    		&my_size,   
    		1,              
    		MPI_INT,
    		size_local.data(),        
    		1,             
    		MPI_INT,
    		mpi_comm);
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::PickLocalGhost(idx_t j, idx_t &ind, Particles* &particle)
{
	if (j < m_Particles_local.Size())
	{
		ind = j;
		particle = &m_Particles_local;
	}
	else
	{
		ind = j - m_Particles_local.Size();
		particle = &m_Particles_ghost;

	}
}

}
