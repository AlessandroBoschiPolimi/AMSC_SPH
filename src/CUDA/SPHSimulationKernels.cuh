#pragma once
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda.h>

#include "Base/SPHSimulation.hpp"
#include "Neighbors/NeighborFinder.cuh"
#include "Neighbors/SpatialHashing.cuh"

#include "Particle.hpp"

namespace cudasph
{
	template <size_t D>
	__global__ void StartKernel()
	{
		// NOTE: this function should be used to create a persistent kernel, that keeps updating without returning to the GPU
		// however the cost to launch a kernel should be negligible, especially when having many particles, so leave this blank for now.

		return;
	}

	template <size_t D>
	__global__ void StepKernel(base::SPHParams params, base::SPHProfiling* dProfiling, ParticlesCuda<D> dParticles)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		dParticles.DSetPositionY(i, dParticles.DPositionY(i) + 0.001);

		return;
	}

	void RunStepKernel(size_t size, base::SPHParams params, base::SPHProfiling* dProfiling, ParticlesCuda<2> dParticles);
	void RunStepKernel(size_t size, base::SPHParams params, base::SPHProfiling* dProfiling, ParticlesCuda<3> dParticles);
}