#pragma once
#ifdef HAS_CUDA

#include "Utility.hpp"
#include "Visualizers/CudaObserver.hpp"
#include "Base/SPHParams.hpp"
#include "Initializers/SimInitializer.hpp"

#include "Neighbors/SpatialHashing.cuh"

#include "SPHSimulationKernels.cuh"

#include <mutex>


namespace cudasph
{

template <size_t D>
struct SpatialHashingStorage;

template <>
struct SpatialHashingStorage<2> { Grid2 Grid; };
template <>
struct SpatialHashingStorage<3> { Grid3 Grid; };


template <size_t D>
class SPHSimulation
{
public:
	using idx_t = Particle<D>::idx_t;
	using vec_t = Particle<D>::vec_t;

	using Particles = ParticleSoA<D>;

	SPHSimulation();
	~SPHSimulation() {}

	void Start();

	Particles& GetParticles(u32 stride);

	void InitializeFluid(const SimInitializer<D, Particles>* init);


	void AddObserver(CudaObserver<D>* obs) {
		obs->Attach(this);
		m_Observers.push_back(obs);
	}


public:
	void SetName(const std::string& name) { m_Name = name; }
	std::string GetName() const { return m_Name; }

	void SetParams(const base::SPHParams& params) { m_Params = params; }
	void SetRestDensity(float val) { m_Params.RestDensity = val; }
	void SetStiffness(float val) { m_Params.Stiffness = val; }
	void SetViscosity(float val) { m_Params.Viscosity = val; }
	void SetTimeStep(float val) { m_Params.TimeStep = val; }
	void SetSmoothingLength(float val) { m_Params.SmoothingLength = val; }
	void SetFinalTime(float val) { m_Params.FinalTime = val; }
	base::SPHParams GetParams() const { return m_Params; }
	float GetRestDensity() const { return m_Params.RestDensity; }
	float GetStiffness() const { return m_Params.Stiffness; }
	float GetViscosity() const { return m_Params.Viscosity; }
	float GetTimeStep() const { return m_Params.TimeStep; }
	float GetSmoothingLength() const { return m_Params.SmoothingLength; }
	float GetFinalTime() const { return m_Params.FinalTime; }

	float GetTime()  const { return m_Time; }
	u64   GetFrame() const { return m_Frame; }
	base::SPHProfiling GetProfiling();

private:
	void Step();

	void NotifyStartFrame() {
		for (auto* o : m_Observers)
			o->OnStartFrame();
	}
	void NotifyEndFrame() {
		for (auto* o : m_Observers)
			o->OnEndFrame();
	}

private:
	base::SPHParams m_Params;

	Particles m_HParticles;
	ParticlesCuda<D> m_DParticles;
	
	bool m_Running = false;
	u64 m_Frame = 0;
	float m_Time = 0.0f;

	base::SPHProfiling m_Profiling;
	//base::SPHProfiling* m_DProfiling;
	std::vector<CudaObserver<D>*> m_Observers;
	std::string m_Name = "Simulation";

	std::mutex m_Mutex;
	size_t m_ParticleCount = 0;

	//NeighborGrid grid;
	SpatialHashingStorage<D> m_Grid;
	GridGPU m_GridGPU;
};


template <size_t D>
inline void SPHSimulation<D>::Start()
{
	std::cout << "Start!\n";

	{
		std::lock_guard<std::mutex> lock(m_Mutex);
		// TODO: send m_Particles to the GPU
		size_t N = m_ParticleCount;
		m_DParticles.Malloc(N);
		m_DParticles.MemcpyFrom(m_HParticles);

		if constexpr (D == 2)
		{
			float2 minBound = make_float2(-0.2, -0.2);
			float2 maxBound = make_float2(1.2, 1.2);
			m_Grid.Grid.init(minBound, maxBound, m_Params.SmoothingLength);
		}
		else
		{
			float3 minBound = make_float3(-0.2, -0.2, -0.2);
			float3 maxBound = make_float3(1.2, 1.2, 1.2);
			m_Grid.Grid.init(minBound, maxBound, m_Params.SmoothingLength);
		}

		mallocGridGPU(m_GridGPU, N, m_Grid.Grid.totalCells);

		//cudaMalloc((void**)&m_DProfiling, sizeof(base::SPHProfiling));
		//cudaMemset(m_DProfiling, 0, sizeof(base::SPHProfiling));

		// initNeighborGrid();

		m_Running = true;
	}

	while (m_Time < m_Params.FinalTime)
		Step();

	{
		std::lock_guard<std::mutex> lock(m_Mutex);
		m_Running = false;
		
		// TODO: free GPU
		m_DParticles.MemcpyTo(m_HParticles);
		m_DParticles.Free();
		freeGridGPU(m_GridGPU);
		//cudaFree(m_DProfiling);
	}
}

template <size_t D>
inline void SPHSimulation<D>::Step()
{
	// NOTE: stepping 1 frame at a time from outside is kinda meaningless, so the UI should just allow pausing and 
	// incrementing FinalTime by a multiple of TimeStep, then resume it; perhaps by adding another param
	// like "Paused" which does std::sleep_for(100ms) in the loop in Start
	// (Pausing the simulation should free the GPU memory)

	NotifyStartFrame();

	cudaError_t cudaStatus;

	RunStepKernels(m_Frame, m_Params, m_DParticles, m_Profiling, m_Grid.Grid, m_GridGPU);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "StepKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		m_Time += m_Params.FinalTime;
	}

	m_Time += m_Params.TimeStep;
	m_Frame++;

	NotifyEndFrame();
}




template <size_t D>
inline SPHSimulation<D>::SPHSimulation()
{
}


template <size_t D>
inline SPHSimulation<D>::Particles& SPHSimulation<D>::GetParticles(u32 stride)
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	std::cout << "CCC\n";
	m_DParticles.MemcpyTo<ParticleSoA<D>>(m_HParticles, stride);
	return m_HParticles;
}
template <size_t D>
inline base::SPHProfiling SPHSimulation<D>::GetProfiling()
{
	//cudaMemcpy((void*)&m_Profiling, m_DProfiling, sizeof(base::SPHProfiling), cudaMemcpyDeviceToHost);
	return m_Profiling;
}

template <size_t D>
inline void SPHSimulation<D>::InitializeFluid(const SimInitializer<D, Particles>* init)
{
	std::cout << "Initializing " << m_Name << '\n';
	std::vector<uptr<Object<2, Particles>>> tmp;
	init->Init(m_HParticles, m_Params.SmoothingLength, tmp);
	m_ParticleCount = m_HParticles.Size();
	std::cout << "Particles: " << m_ParticleCount << '\n';
}

}

#endif