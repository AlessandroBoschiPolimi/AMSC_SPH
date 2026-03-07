#include "Neighbors/NeighborFinder.cuh"
#include "Base/SPHSimulation.hpp"

#include "Neighbors/SpatialHashing.cuh"

#include "Utility.hpp"
#include <mutex>

#include "SPHSimulationKernels.cuh"


namespace cudasph
{

template <size_t D, ParticleSet<D> Particles = ParticleSoA<D>>
class SPHSimulation : public base::SPHSimulation<D, Particles>
{
public:
	using idx_t = Particle<D>::idx_t;
	using vec_t = Particle<D>::vec_t;


	SPHSimulation();
	virtual ~SPHSimulation() override {}

	void Start() override;


	const Particles& GetParticles() const override;
	Particles& GetParticles() override;
	const std::vector<std::vector<idx_t>>& GetNeighbors() const override;
	std::vector<std::vector<idx_t>>& GetNeighbors() override;

	void InitializeFluid(const SimInitializer<D, Particles>* init) override;


protected:
	void Step() override;

private:
	mutable Particles m_HParticles;
	ParticlesCuda<D> m_DParticles;

	base::SPHProfiling* m_DProfiling;

	std::vector<std::vector<idx_t>> m_HNeighbors_UNUSED;

	bool m_Running = false;

	mutable std::mutex m_Mutex;
};


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

		cudaMalloc((void**)&m_DProfiling, sizeof(base::SPHProfiling));
		cudaMemset(m_DProfiling, 0, sizeof(base::SPHProfiling));

		// initNeighborGrid();

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
		cudaFree(m_DProfiling);
	}
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::Step()
{
	// NOTE: stepping 1 frame at a time from outside is kinda meaningless, so the UI should just allow pausing and 
	// incrementing FinalTime by a multiple of TimeStep, then resume it; perhaps by adding another param
	// like "Paused" which does std::sleep_for(100ms) in the loop in Start
	// (Pausing the simulation should free the GPU memory)

	this->NotifyStartFrame();

	cudaError_t cudaStatus;

	RunStepKernel(m_HParticles.Size(), this->m_Params, m_DProfiling, m_DParticles);

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
inline SPHSimulation<D, Particles>::SPHSimulation()
	: base::SPHSimulation<D, Particles>()
{
}


template <size_t D, ParticleSet<D> Particles>
inline Particles& SPHSimulation<D, Particles>::GetParticles()
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	m_DParticles.MemcpyTo<ParticleSoA<D>>(m_HParticles);
	cudaMemcpy((void*)&this->m_Profiling, m_DProfiling, sizeof(base::SPHProfiling), cudaMemcpyDeviceToHost);
	return m_HParticles;
}
template <size_t D, ParticleSet<D> Particles>
inline const Particles& SPHSimulation<D, Particles>::GetParticles() const
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	m_DParticles.MemcpyTo<ParticleSoA<D>>(m_HParticles);
	cudaMemcpy((void*)&this->m_Profiling, m_DProfiling, sizeof(base::SPHProfiling), cudaMemcpyDeviceToHost);
	return m_HParticles;
}
template <size_t D, ParticleSet<D> Particles>
inline std::vector<std::vector<typename SPHSimulation<D, Particles>::idx_t>>& SPHSimulation<D, Particles>::GetNeighbors()
{
	return m_HNeighbors_UNUSED;
}
template <size_t D, ParticleSet<D> Particles>
inline const std::vector<std::vector<typename SPHSimulation<D, Particles>::idx_t>>& SPHSimulation<D, Particles>::GetNeighbors() const
{
	return m_HNeighbors_UNUSED;
}

template <size_t D, ParticleSet<D> Particles>
inline void SPHSimulation<D, Particles>::InitializeFluid(const SimInitializer<D, Particles>* init)
{
	std::cout << "Initializing " << this->m_Name << '\n';
	init->Init(m_HParticles, this->m_Params.SmoothingLength, this->m_Objects);
	std::cout << "Particles: " << m_HParticles.Size() << '\n';
}

}