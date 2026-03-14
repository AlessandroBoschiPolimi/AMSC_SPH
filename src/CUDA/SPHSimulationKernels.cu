#ifdef HAS_CUDA
#include "SPHSimulationKernels.cuh"

namespace cudasph
{
	__device__ ParticlesCuda<2>* DParticles2;
	__device__ ParticlesCuda<3>* DParticles3;
	__device__ base::SPHParams DParams;
	__device__ float DSL2;
	__device__ Grid2 DGrid2;
	__device__ Grid3 DGrid3;
	__device__ GridGPU DGridGPU;

	__device__ float KernelValue3(const typename Particle<3>::vec_t& r)
	{
		float h = DParams.SmoothingLength;

		float r2 = r.x * r.x + r.y * r.y + r.z * r.z;

		float q = sqrt(r2) / h;
		float out = 0.0f;
		if (q < 1.0f)
		{
			float x1 = 2.0f - q;
			float x2 = 1.0f - q;
			out = x1 * x1 * x1 - 4.0f * x2 * x2 * x2;
		}
		else if (q >= 1.0f && q < 2.0f)
		{
			float x1 = 2.0f - q;
			out = x1 * x1 * x1;
		}

		float alpha = 1.0f / (4.0f * 3.14 * h * h * h);

		return alpha * out;
	}
	__device__ Particle<3>::vec_t KernelGradient3(const typename Particle<3>::vec_t& r)
	{
		float h = DParams.SmoothingLength;

		float r2 = r.x * r.x + r.y * r.y + r.z * r.z;

		float q = sqrt(r2) / h;
		float dW_dq = 0.0f;
		if (q < 1.0f)
		{
			float x1 = 2.0f - q;
			float x2 = 1.0f - q;
			dW_dq = -3.0f * x1 * x1 + 12.0f * x2 * x2;
		}
		else if (q >= 1.0f && q < 2.0f)
		{
			float x1 = 2.0f - q;
			dW_dq = -3.0f * x1 * x1;
		}

		float alpha = 1.0f / (4.0f * 3.14 * h * h * h);

		float prefac = alpha * dW_dq / (sqrt(r2) * h);
		return prefac * r;
	}

	struct BoundaryPsiCallback3
	{
		using vec_t = Particle<3>::vec_t;
		float V = 0;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			if (DParticles3->Types[j] == SOLID)
				V += KernelValue3(r);
		}
	};
	struct DensityCallback3
	{
		using vec_t = Particle<3>::vec_t;
		float density;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			float W_ij = KernelValue3(r);
			if (DParticles3->Types[j] == FLUID)
				density += DParticles3->Masss[j] * W_ij;
			else // Boundary handling
				density += DParticles3->BoundaryPsis[j] * W_ij;
		}
	};
	struct AccelerationViscosityCallback3
	{
		using vec_t = Particle<3>::vec_t;

		vec_t pi;
		vec_t vi;
		float di;
		vec_t acc;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			vec_t DW_ij = KernelGradient3(r);
			vec_t v_ij = vi - DParticles3->DVelocity(j);
			vec_t x_ij = pi - DParticles3->DPosition(j);
			float prod = Dot(x_ij, v_ij);
			float numerator = ((prod > 0) ? prod : 0);

			if (DParticles3->Types[j] == FLUID)
			{
				numerator *= DParticles3->Masss[j];
				float denominator = (di + DParticles3->Densitys[j]) * (Dot(x_ij, x_ij) + 0.01 * DSL2);
				acc += (2 * DParams.Viscosity * (numerator / denominator)) * DW_ij;
			}
			else
			{
				float denominator = 2 * di * (Dot(x_ij, x_ij) + 0.01 * DSL2);
				acc += (DParams.ViscosityRigid * (numerator / denominator)) * DParticles3->BoundaryPsis[j] * DW_ij;
			}
		}
	};
	struct AccelerationPressureCallback3
	{
		using vec_t = Particle<3>::vec_t;

		vec_t pi;
		float pri;
		float ddi;

		vec_t acc;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			float factor = (DParticles3->Types[j] == SOLID) ? (DParticles3->BoundaryPsis[j] / 2.0f) : DParticles3->Masss[j];
			u32 z = (DParticles3->Types[j] == SOLID) ? i : j;
			vec_t DW_ij = KernelGradient3(r);
			auto dz = DParticles3->Densitys[z];

			acc -= (pri / ddi + DParticles3->Pressures[z] / (dz * dz)) * factor * DW_ij;
		}
	};


	__device__ void DComputeBoundaryPsi3(u32 i)
	{
		/*
		 * Computes 'mass' of the boundary particles used
		 * to implement collisions
		 */

		BoundaryPsiCallback3 callback;

		
		float3 pos_i = make_float3(DParticles3->DPositionX(i), DParticles3->DPositionY(i), DParticles3->DPositionZ(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles3->Xs, DParticles3->Ys, DParticles3->Zs,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid3,
			callback
		);

		// Clamp the values in case the volume is too small
		DParticles3->DSetBoundaryPsi(i, (callback.V > 1.0f) ? DParams.RestDensity / callback.V : 0);
	}
	__device__ void DComputeDensity3(u32 i)
	{
		/*
		 * Simple function to compute density
		 * Takes initial value of density produced by itself
		 * For security, clamps the value in the end to avoid disappearing particles
		 */

		DensityCallback3 callback;

		callback.density = DParticles3->DMass(i) * KernelValue3({ 0, 0, 0 });
		float3 pos_i = make_float3(DParticles3->DPositionX(i), DParticles3->DPositionY(i), DParticles3->DPositionZ(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles3->Xs, DParticles3->Ys, DParticles3->Zs,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid3,
			callback
		);

		DParticles3->DSetDensity(i, callback.density);
	}
	__device__ void DComputePressure3(u32 i)
	{
		/*
		 * (Andrew) Tait equation
		 * Stiffness constant is user defined
		 */

		float pressure = max(
			DParams.Stiffness * (pow(DParticles3->DDensity(i) / (DParams.RestDensity), 7.0f) - 1),
			0.0f
		);
		DParticles3->DSetPressure(i, pressure);
	}
	__device__ void DComputeAccelerationViscosity3(u32 i)
	{
		/*
		 * Computes acceleration dues to viscosity from [1]
		 * [2] has typos in that formula
		 */

		using vec_t = Particle<3>::vec_t;


		AccelerationViscosityCallback3 callback;
		callback.pi = DParticles3->DPosition(i);
		callback.vi = DParticles3->DVelocity(i);
		callback.di = DParticles3->DDensity(i);

		callback.acc = vec_t{ 0, 0, 0 };
		float3 pos_i = make_float3(DParticles3->DPositionX(i), DParticles3->DPositionY(i), DParticles3->DPositionZ(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles3->Xs, DParticles3->Ys, DParticles3->Zs,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid3,
			callback
		);

		DParticles3->DSetAccelerationViscosity(i, callback.acc);
	}
	__device__ void DComputeAccelerationPressure3(u32 i)
	{
		/*
		 * Computes acceleration dues to viscosity from [2]
		 * Boundary handling according to [2]
		 */

		using vec_t = Particle<3>::vec_t;
		using idx_t = Particle<3>::idx_t;

		AccelerationPressureCallback3 callback;
		
		callback.pi = DParticles3->DPosition(i);
		callback.pri = DParticles3->DPressure(i);
		auto di = DParticles3->DDensity(i);
		callback.ddi = di * di;

		callback.acc = vec_t{ 0, 0, 0 };
		float3 pos_i = make_float3(DParticles3->DPositionX(i), DParticles3->DPositionY(i), DParticles3->DPositionZ(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles3->Xs, DParticles3->Ys, DParticles3->Zs,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid3,
			callback
		);

		DParticles3->DSetAccelerationPressure(i, callback.acc);
	}
	__device__ inline void DUpdatePositionInitial3(u32 i)
	{
		// Update position due to viscosity and gravity
		auto pi = DParticles3->DPosition(i);
		DParticles3->DSetPosition(i, pi + DParams.TimeStep * DParticles3->DVelocity(i));
	}
	__device__ inline void DUpdatePositionIteration3(u32 i)
	{
		// Update position due to pressure
		auto pi = DParticles3->DPosition(i);
		pi += DParams.TimeStep * DParams.TimeStep *
			DParticles3->DAccelerationPressure(i);
		DParticles3->DSetPosition(i, pi);
	}
	__device__ inline void DUpdateVelocityInitial3(u32 i)
	{
		// Update velocity due to viscosity and gravity
		auto vi = DParticles3->DVelocity(i);
		vi += DParams.TimeStep * (DParticles3->DAccelerationGravity(i) + DParticles3->DAccelerationViscosity(i));
		DParticles3->DSetVelocity(i, vi);
	}
	__device__ inline void DUpdateVelocityIteration3(u32 i)
	{
		// Update velocity due to pressure
		auto vi = DParticles3->DVelocity(i);
		vi += DParams.TimeStep * DParticles3->DAccelerationPressure(i);
		DParticles3->DSetVelocity(i, vi);
	}


	__global__ void ComputeDensityOrBoundaryPsiKernel3()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles3->Count) return;

		if (DParticles3->DType(i) == FLUID)
			DComputeDensity3(i);
		else
			DComputeBoundaryPsi3(i);
	}
	__global__ void ComputeAccelerationViscosityKernel3()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles3->Count) return;

		if (DParticles3->DType(i) == FLUID)
			DComputeAccelerationViscosity3(i);
	}
	__global__ void UpdateVelocityPositionInitialKernel3()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles3->Count) return;

		if (DParticles3->DType(i) == FLUID)
		{
			DUpdateVelocityInitial3(i);
			DUpdatePositionInitial3(i);
		}
	}

	__global__ void ComputeDensityAndPressureKernel3()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles3->Count) return;

		if (DParticles3->DType(i) == FLUID)
		{
			DComputeDensity3(i);
			DComputePressure3(i);
		}
	}
	__global__ void ComputeAccelerationPressureKernel3()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles3->Count) return;

		if (DParticles3->DType(i) == FLUID)
			DComputeAccelerationPressure3(i);
	}
	__global__ void UpdateVelocityPositionIterationKernel3()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles3->Count) return;

		if (DParticles3->DType(i) == FLUID)
		{
			DUpdatePositionIteration3(i);
			DUpdateVelocityIteration3(i);
		}
	}






	__device__ float KernelValue2(const typename Particle<2>::vec_t& r)
	{
		float h = DParams.SmoothingLength;

		float r2 = r.x * r.x + r.y * r.y;

		float q = sqrt(r2) / h;
		float out = 0.0f;
		if (q < 1.0f)
		{
			float x1 = 2.0f - q;
			float x2 = 1.0f - q;
			out = x1 * x1 * x1 - 4.0f * x2 * x2 * x2;
		}
		else if (q >= 1.0f && q < 2.0f)
		{
			float x1 = 2.0f - q;
			out = x1 * x1 * x1;
		}

		float alpha = 5.0f / (14.0f * 3.14 * h * h);

		return alpha * out;
	}
	__device__ Particle<2>::vec_t KernelGradient2(const typename Particle<2>::vec_t& r)
	{
		float h = DParams.SmoothingLength;

		float r2 = r.x * r.x + r.y * r.y;

		float q = sqrt(r2) / h;
		float dW_dq = 0.0f;
		if (q < 1.0f)
		{
			float x1 = 2.0f - q;
			float x2 = 1.0f - q;
			dW_dq = -3.0f * x1 * x1 + 12.0f * x2 * x2;
		}
		else if (q >= 1.0f && q < 2.0f)
		{
			float x1 = 2.0f - q;
			dW_dq = -3.0f * x1 * x1;
		}

		float alpha = 5.0f / (14.0f * 3.14 * h * h);

		float prefac = alpha * dW_dq / (sqrt(r2) * h);
		return prefac * r;
	}

	struct BoundaryPsiCallback2
	{
		using vec_t = Particle<2>::vec_t;
		float V = 0;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			if (DParticles2->Types[j] == SOLID)
				V += KernelValue2(r);
		}
	};
	struct DensityCallback2
	{
		using vec_t = Particle<2>::vec_t;
		float density;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			float W_ij = KernelValue2(r);
			if (DParticles2->Types[j] == FLUID)
				density += DParticles2->Masss[j] * W_ij;
			else // Boundary handling
				density += DParticles2->BoundaryPsis[j] * W_ij;
		}
	};
	struct AccelerationViscosityCallback2
	{
		using vec_t = Particle<2>::vec_t;

		vec_t pi;
		vec_t vi;
		float di;
		vec_t acc;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			vec_t DW_ij = KernelGradient2(r);
			vec_t v_ij = vi - DParticles2->DVelocity(j);
			vec_t x_ij = pi - DParticles2->DPosition(j);
			float prod = Dot(x_ij, v_ij);
			float numerator = ((prod > 0) ? prod : 0);

			if (DParticles2->Types[j] == FLUID)
			{
				numerator *= DParticles2->Masss[j];
				float denominator = (di + DParticles2->Densitys[j]) * (Dot(x_ij, x_ij) + 0.01 * DSL2);
				acc += (2 * DParams.Viscosity * (numerator / denominator)) * DW_ij;
			}
			else
			{
				float denominator = 2 * di * (Dot(x_ij, x_ij) + 0.01 * DSL2);
				acc += (DParams.ViscosityRigid * (numerator / denominator)) * DParticles2->BoundaryPsis[j] * DW_ij;
			}
		}
	};
	struct AccelerationPressureCallback2
	{
		using vec_t = Particle<2>::vec_t;

		vec_t pi;
		float pri;
		float ddi;

		vec_t acc;

		__device__ void operator()(int i, int j, const vec_t& r, float r2)
		{
			float factor = (DParticles2->Types[j] == SOLID) ? (DParticles2->BoundaryPsis[j] / 2.0f) : DParticles2->Masss[j];
			u32 z = (DParticles2->Types[j] == SOLID) ? i : j;
			vec_t DW_ij = KernelGradient2(r);
			auto dz = DParticles2->Densitys[z];

			acc -= (pri / ddi + DParticles2->Pressures[z] / (dz * dz)) * factor * DW_ij;
		}
	};


	__device__ void DComputeBoundaryPsi2(u32 i)
	{
		/*
		 * Computes 'mass' of the boundary particles used
		 * to implement collisions
		 */

		BoundaryPsiCallback2 callback;

		float2 pos_i = make_float2(DParticles2->DPositionX(i), DParticles2->DPositionY(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles2->Xs, DParticles2->Ys,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid2,
			callback
		);

		// Clamp the values in case the volume is too small
		DParticles2->DSetBoundaryPsi(i, (callback.V > 1.0f) ? DParams.RestDensity / callback.V : 0);
	}
	__device__ void DComputeDensity2(u32 i)
	{
		/*
		 * Simple function to compute density
		 * Takes initial value of density produced by itself
		 * For security, clamps the value in the end to avoid disappearing particles
		 */

		DensityCallback2 callback;

		callback.density = DParticles2->DMass(i) * KernelValue2({ 0, 0 });
		float2 pos_i = make_float2(DParticles2->DPositionX(i), DParticles2->DPositionY(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles2->Xs, DParticles2->Ys,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid2,
			callback
		);

		DParticles2->DSetDensity(i, callback.density);
	}
	__device__ void DComputePressure2(u32 i)
	{
		/*
		 * (Andrew) Tait equation
		 * Stiffness constant is user defined
		 */

		float pressure = max(
			DParams.Stiffness * (pow(DParticles2->DDensity(i) / (DParams.RestDensity), 7.0f) - 1),
			0.0f
		);
		DParticles2->DSetPressure(i, pressure);
	}
	__device__ void DComputeAccelerationViscosity2(u32 i)
	{
		/*
		 * Computes acceleration dues to viscosity from [1]
		 * [2] has typos in that formula
		 */

		using vec_t = Particle<2>::vec_t;


		AccelerationViscosityCallback2 callback;
		callback.pi = DParticles2->DPosition(i);
		callback.vi = DParticles2->DVelocity(i);
		callback.di = DParticles2->DDensity(i);

		callback.acc = vec_t{ 0, 0 };
		float2 pos_i = make_float2(DParticles2->DPositionX(i), DParticles2->DPositionY(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles2->Xs, DParticles2->Ys,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid2,
			callback
		);

		DParticles2->DSetAccelerationViscosity(i, callback.acc);
	}
	__device__ void DComputeAccelerationPressure2(u32 i)
	{
		/*
		 * Computes acceleration dues to viscosity from [2]
		 * Boundary handling according to [2]
		 */

		using vec_t = Particle<2>::vec_t;
		using idx_t = Particle<2>::idx_t;

		AccelerationPressureCallback2 callback;

		callback.pi = DParticles2->DPosition(i);
		callback.pri = DParticles2->DPressure(i);
		auto di = DParticles2->DDensity(i);
		callback.ddi = di * di;

		callback.acc = vec_t{ 0, 0 };
		float2 pos_i = make_float2(DParticles2->DPositionX(i), DParticles2->DPositionY(i));
		forEachNeighbor(
			i,
			pos_i,
			DSL2,
			DParticles2->Xs, DParticles2->Ys,
			DGridGPU.cellStart,
			DGridGPU.cellEnd,
			DGrid2,
			callback
		);

		DParticles2->DSetAccelerationPressure(i, callback.acc);
	}
	__device__ inline void DUpdatePositionInitial2(u32 i)
	{
		// Update position due to viscosity and gravity
		auto pi = DParticles2->DPosition(i);
		DParticles2->DSetPosition(i, pi + DParams.TimeStep * DParticles2->DVelocity(i));
	}
	__device__ inline void DUpdatePositionIteration2(u32 i)
	{
		// Update position due to pressure
		auto pi = DParticles2->DPosition(i);
		pi += DParams.TimeStep * DParams.TimeStep *
			DParticles2->DAccelerationPressure(i);
		DParticles2->DSetPosition(i, pi);
	}
	__device__ inline void DUpdateVelocityInitial2(u32 i)
	{
		// Update velocity due to viscosity and gravity
		auto vi = DParticles2->DVelocity(i);
		vi += DParams.TimeStep * (DParticles2->DAccelerationGravity(i) + DParticles2->DAccelerationViscosity(i));
		DParticles2->DSetVelocity(i, vi);
	}
	__device__ inline void DUpdateVelocityIteration2(u32 i)
	{
		// Update velocity due to pressure
		auto vi = DParticles2->DVelocity(i);
		vi += DParams.TimeStep * DParticles2->DAccelerationPressure(i);
		DParticles2->DSetVelocity(i, vi);
	}


	__global__ void ComputeDensityOrBoundaryPsiKernel2()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles2->Count) return;

		if (DParticles2->DType(i) == FLUID)
			DComputeDensity2(i);
		else
			DComputeBoundaryPsi2(i);
	}
	__global__ void ComputeAccelerationViscosityKernel2()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles2->Count) return;

		if (DParticles2->DType(i) == FLUID)
			DComputeAccelerationViscosity2(i);
	}
	__global__ void UpdateVelocityPositionInitialKernel2()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles2->Count) return;

		if (DParticles2->DType(i) == FLUID)
		{
			DUpdateVelocityInitial2(i);
			DUpdatePositionInitial2(i);
		}
	}

	__global__ void ComputeDensityAndPressureKernel2()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles2->Count) return;

		if (DParticles2->DType(i) == FLUID)
		{
			DComputeDensity2(i);
			DComputePressure2(i);
		}
	}
	__global__ void ComputeAccelerationPressureKernel2()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles2->Count) return;

		if (DParticles2->DType(i) == FLUID)
			DComputeAccelerationPressure2(i);
	}
	__global__ void UpdateVelocityPositionIterationKernel2()
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= DParticles2->Count) return;

		if (DParticles2->DType(i) == FLUID)
		{
			DUpdatePositionIteration2(i);
			DUpdateVelocityIteration2(i);
		}
	}




	void RunStepKernels(
		size_t frame, base::SPHParams params, 
		ParticlesCuda<2> dParticles, base::SPHProfiling& profiling,
		Grid2& grid, GridGPU& gridGPU)
	{
		int threads = 256;
		int blocks = (dParticles.Size() + threads - 1) / threads;
		size_t sharedMem = sizeof(float2) * threads;

		float sl2 = params.SmoothingLength * params.SmoothingLength;

		std::cout << "BANANA 1\n";
		
		cudaError_t cudaStatus;

		if (frame == 0)
		{
			ParticlesCuda<2>* devPtr;
			cudaMalloc(&devPtr, sizeof(ParticlesCuda<2>));
			cudaMemcpy(devPtr, &dParticles, sizeof(dParticles), cudaMemcpyHostToDevice);

			cudaMemcpyToSymbol(DParticles2, &devPtr, sizeof(devPtr));
		}
		cudaMemcpyToSymbol(DParams, &params, sizeof(base::SPHParams));
		cudaMemcpyToSymbol(DSL2, &sl2, sizeof(float));
		cudaMemcpyToSymbol(DGrid2, &grid, sizeof(Grid2));
		cudaMemcpyToSymbol(DGridGPU, &gridGPU, sizeof(GridGPU));
		
		cudaDeviceSynchronize();

		{
			stdc::time_point<stdclock> start;
			start = stdclock::now();

			std::cout << "BANANA 2\n";
			buildGrid(dParticles, grid, gridGPU);

			cudaError_t err = cudaDeviceSynchronize();
			if (err != cudaSuccess)
				printf("%s\n", cudaGetErrorString(err));

			profiling.Neighbors = stdclock::now() - start;
		}
		{
			stdc::time_point<stdclock> start;
			start = stdclock::now();

			std::cout << "BANANA 3\n";
			if (frame == 0)
			{
				ComputeDensityOrBoundaryPsiKernel2<<<blocks, threads, sharedMem>>>();

				cudaError_t err = cudaDeviceSynchronize();
				if (err != cudaSuccess)
					printf("%s\n", cudaGetErrorString(err));

				std::cout << "BANANA 31\n";

				ComputeAccelerationViscosityKernel2<<<blocks, threads, sharedMem>>>();

				err = cudaDeviceSynchronize();
				if (err != cudaSuccess)
					printf("%s\n", cudaGetErrorString(err));
			}

			std::cout << "BANANA 32\n";

			UpdateVelocityPositionInitialKernel2<<<blocks, threads, sharedMem>>>();

			cudaError_t err = cudaDeviceSynchronize();
			if (err != cudaSuccess)
				printf("%s\n", cudaGetErrorString(err));
		
			profiling.Initialize = stdclock::now() - start;
		}

		{
			stdc::time_point<stdclock> start;
			start = stdclock::now();

			std::cout << "BANANA 4\n";
			ComputeDensityAndPressureKernel2<<<blocks, threads, sharedMem>>>();
			ComputeAccelerationPressureKernel2<<<blocks, threads, sharedMem>>>();
			UpdateVelocityPositionIterationKernel2<<<blocks, threads, sharedMem>>>();

			cudaError_t err = cudaDeviceSynchronize();
			if (err != cudaSuccess)
				printf("%s\n", cudaGetErrorString(err));

			profiling.IterativePressure = stdclock::now() - start;
		}

		std::cout << "BANANA 5\n";
	}
	void RunStepKernels(
		size_t frame, base::SPHParams params, 
		ParticlesCuda<3> dParticles, base::SPHProfiling& profiling,
		Grid3& grid, GridGPU& gridGPU)
	{
		int threads = 256;
		int blocks = (dParticles.Size() + threads - 1) / threads;
		size_t sharedMem = sizeof(float3) * threads;

		float sl2 = params.SmoothingLength * params.SmoothingLength;

		if (frame == 0)
		{
			ParticlesCuda<3>* devPtr;
			cudaMalloc(&devPtr, sizeof(ParticlesCuda<3>));
			cudaMemcpy(devPtr, &dParticles, sizeof(dParticles), cudaMemcpyHostToDevice);

			cudaMemcpyToSymbol(DParticles3, &devPtr, sizeof(devPtr));
		}
		cudaMemcpyToSymbol(DParams, &params, sizeof(base::SPHParams));
		cudaMemcpyToSymbol(DSL2, &sl2, sizeof(float));
		cudaMemcpyToSymbol(DGrid3, &grid, sizeof(Grid3));
		cudaMemcpyToSymbol(DGridGPU, &gridGPU, sizeof(GridGPU));

		{
			stdc::time_point<stdclock> start;
			start = stdclock::now();
			
			buildGrid(dParticles, grid, gridGPU);

			profiling.Neighbors = stdclock::now() - start;
		}
		{
			stdc::time_point<stdclock> start;
			start = stdclock::now();

			if (frame == 0)
			{
				ComputeDensityOrBoundaryPsiKernel3<<<blocks, threads, sharedMem>>>();
				ComputeAccelerationViscosityKernel3<<<blocks, threads, sharedMem>>>();
			}

			UpdateVelocityPositionInitialKernel3<<<blocks, threads, sharedMem>>>();

			profiling.Initialize = stdclock::now() - start;
		}

		{
			stdc::time_point<stdclock> start;
			start = stdclock::now();

			ComputeDensityAndPressureKernel3<<<blocks, threads, sharedMem>>>();
			ComputeAccelerationPressureKernel3<<<blocks, threads, sharedMem>>>();
			UpdateVelocityPositionIterationKernel3<<<blocks, threads, sharedMem>>>();

			profiling.IterativePressure = stdclock::now() - start;
		}
	}
}
#endif