#pragma once
#ifdef HAS_CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda.h>

#include "Base/SPHParams.hpp"
#include "Neighbors/SpatialHashing.cuh"

#include "Particle.hpp"

template <size_t D>
struct ParticleCuda;

namespace cudasph
{
	template <size_t D>
	__global__ void StartKernel()
	{
		// NOTE: this function should be used to create a persistent kernel, that keeps updating without returning to the GPU
		// however the cost to launch a kernel should be negligible, especially when having many particles, so leave this blank for now.

		return;
	}


	void RunStepKernels(size_t frame, base::SPHParams params, ParticlesCuda<2> dParticles, base::SPHProfiling& profiling, Grid2& grid, GridGPU& gridGPU);
	void RunStepKernels(size_t frame, base::SPHParams params, ParticlesCuda<3> dParticles, base::SPHProfiling& profiling, Grid3& grid, GridGPU& gridGPU);
}
#endif