#ifdef HAS_CUDA
#include "SPHSimulationKernels.cuh"

namespace cudasph
{
	void RunStepKernel(size_t size, base::SPHParams params, base::SPHProfiling* dProfiling, ParticlesCuda<2> dParticles)
	{
		int threads = 256;
		int blocks = (size + threads - 1) / threads;

		StepKernel<<<blocks, threads>>>(params, dProfiling, dParticles);
	}
	void RunStepKernel(size_t size, base::SPHParams params, base::SPHProfiling* dProfiling, ParticlesCuda<3> dParticles)
	{
		int threads = 256;
		int blocks = (size + threads - 1) / threads;

		StepKernel<<<blocks, threads>>>(params, dProfiling, dParticles);
	}
}
#endif