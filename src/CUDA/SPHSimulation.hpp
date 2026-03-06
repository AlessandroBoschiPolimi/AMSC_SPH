#include "Neighbors/NeighborFinder.hpp"
#include "Base/SPHSimulation.hpp"


namespace cuda
{

template <size_t D>
__global__ void StartKernel();
template <size_t D>
__global__ void StepKernel();

template <size_t D, ParticleSet<D> Particles = ParticleSoA<D>>
class SPHSimulation : public base::SPHSimulation<D, Particles>
{
public:
	using idx_t = Particle<D>::idx_t;
	using vec_t = Particle<D>::vec_t;

	friend class NeighborFinder<D, Particles>;


	SPHSimulation(NeighborFinder<D, Particles>* nf);
	virtual ~SPHSimulation() override {}

	void Start() override;


	const Particles& GetParticles() const override;
	Particles& GetParticles() override;
	const std::vector<std::vector<idx_t>>& GetNeighbors() const override;
	std::vector<std::vector<idx_t>>& GetNeighbors() override;

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
	Particles m_HParticles;
	ParticlesCuda<D> m_DParticles;
	std::vector<std::vector<idx_t>> m_HNeighbors;

	NeighborFinder<D, Particles>* m_NeighborFinder = nullptr;

	Kernel<D> W_Ker;

	bool m_Running = false;
	std::mutex m_Lock;
};


template <size_t D, ParticleSet<D> Particles>
inline SPHSimulation<D, Particles>::SPHSimulation(NeighborFinder<D, Particles>* nf)
	: base::SPHSimulation<D, Particles>(), m_NeighborFinder(nf), W_Ker(this->m_Params.SmoothingLength / 2.0f)
{
}


template <size_t D, ParticleSet<D> Particles>
inline Particles& SPHSimulation<D, Particles>::GetParticles()
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	if (m_Running)
		cudaMemcpy(m_HParticles.data(), m_DParticles, N * sizeof(float), cudaMemcpyDeviceToHost);
	return m_HParticles;
}
template <size_t D, ParticleSet<D> Particles>
inline const Particles& SPHSimulation<D, Particles>::GetParticles() const
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	if (m_Running)
		cudaMemcpy(m_HParticles.data(), m_DParticles, N * sizeof(float), cudaMemcpyDeviceToHost);
	return m_HParticles;
}
template <size_t D, ParticleSet<D> Particles>
inline std::vector<std::vector<typename SPHSimulation<D, Particles>::idx_t>>& SPHSimulation<D, Particles>::GetNeighbors()
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	if (m_Running)
		cudaMemcpy(m_HNeighbors.data(), m_DNeighbors, N * sizeof(float), cudaMemcpyDeviceToHost);
	return m_HNeighbors;
}
template <size_t D, ParticleSet<D> Particles>
inline const std::vector<std::vector<typename SPHSimulation<D, Particles>::idx_t>>& SPHSimulation<D, Particles>::GetNeighbors() const
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	if (m_Running)
		cudaMemcpy(m_HNeighbors.data(), m_DNeighbors, N * sizeof(float), cudaMemcpyDeviceToHost);
	return m_HNeighbors;
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Start()
{
	std::cout << "Start!\n";

	{
		std::lock_guard<std::mutex> lock(m_Mutex);
		// TODO: send m_Particles to the GPU
		size_t N = m_HParticles.Size();
		m_DParticles.Malloc(N);
		m_DParticles.MemcpyFrom(m_HParticles);
		cudaMalloc((void**)&m_DNeighbors, N * sizeof(float));

		cudaMalloc(&grid.d_cellHash, N * sizeof(uint32_t));
		cudaMalloc(&grid.d_particleIndex, N * sizeof(uint32_t));

		size_t numCells = N * 2;
		cudaMalloc(&grid.d_cellStart, numCells * sizeof(uint32_t));
		cudaMalloc(&grid.d_cellEnd, numCells * sizeof(uint32_t));

		cudaMalloc(&grid.d_sortedPos, N * sizeof(float3));

		m_Running = true;
	}

	while (this->m_Time < this->m_Params.FinalTime)
		Step();

	{
		std::lock_guard<std::mutex> lock(m_Mutex);
		m_Running = false;
		// TODO: free GPU
		m_DParticles.MemcpyTo(m_HParticles);
		m_DParticles.Free();
		cudaFree(m_DNeighbors);
	}
}
template <size_t D, ParticleSet<D> Particles>
inline __global__ void StartKernel()
{
	// NOTE: this function should be used to create a persistent kernel, that keeps updating without returning to the GPU
	// however the cost to launch a kernel should be negligible, especially when having many particles, so leave this blank for now.

	return;
}
template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Step()
{
	// NOTE: stepping 1 frame at a time from outside is kinda meaningless, so the UI should just allow pausing and 
	// incrementing FinalTime by a multiple of TimeStep, then resume it; perhaps by adding another param
	// like "Paused" which does std::sleep_for(100ms) in the loop in Start
	// (Pausing the simulation should free the GPU memory)
	// TODO: make this function protected.

	this->NotifyStartFrame();

	cudaError_t cudaStatus;

	// TODO: to support adding/removing particles we have to remove on host, allocate enough memory and then add on device
	// to remove particles do the following, or add another array of bools to mark them for removal at .
	/* remove particles
	thrust::device_ptr<float> begin(d_data);
	thrust::device_ptr<float> end = begin + size;

	end = thrust::remove_if(begin, end, [] __device__(float x) {
		return x < 0.0f;  // remove negative elements
	});

	size = end - begin;

	if (neighbors_size != particles_size)
	{
		// reallocate
	}
	*/

	StepKernel<<<>>>(this->m_Params, m_DProfiling, m_DParticles, m_DNeighbors);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "StepKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		this->m_Time += this->m_Params.FinalTime;
	}

	this->m_Time += this->m_Params.TimeStep;
	this->m_Frame++;

	this->NotifyEndFrame();

}
template <size_t D, ParticleSet<D> Particles>
inline __global__ void StepKernel(base::SPHParams params, base::SPHProfiling dProfiling, ParticlesCuda* dParticles, Particle<D>::idx_t** dNeighbors)
{
	/* add particles
	int idx = threadIdx.x + blockIdx.x * blockDim.x;

	if (should_add(idx)) {
		int pos = atomicAdd(&d_size, 1); // or = idx if position marked for removal.
		data[pos] = compute_value(idx);
	}
	*/

	//for (auto& obj : m_Objects)
	//	obj->OnFrameStart();

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
	return;
}


template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::FindAllNeighbors()
{
	BuildSpatialGrid(grid, d_pos);

	NeighborSearchKernel << <blocks, threads >> > (
		N,
		grid.d_sortedPos,
		smoothingLength,
		smoothingLength * smoothingLength,

		grid.d_particleIndex,
		grid.d_cellStart,
		grid.d_cellEnd,

		d_neighborList,
		d_neighborCount,
		MAX_NEIGHBORS
		);

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
	vec_t acc = vec_t{ 0, 0 };

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

	vec_t acc = vec_t{ 0, 0 };

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

	for (auto& obj : this->m_Objects)
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