#pragma once
#ifdef HAS_CUDA
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <device_launch_parameters.h>
#include <cuda.h>

__global__ void SampleStrideKernel(
	size_t sample_size, size_t stride,
	const float* Xs, const float* Ys, const float* Zs,
	const float* VXs, const float* VYs, const float* VZs,
	const float* AX_gravs, const float* AY_gravs, const float* AZ_gravs,
	const float* AX_presss, const float* AY_presss, const float* AZ_presss,
	const float* AX_viscs, const float* AY_viscs, const float* AZ_viscs,
	const float* Masss, const float* Densitys, const float* Pressures,
	const uint8_t* Types, const float* BoundaryPsis,
	float* Xs_out, float* Ys_out, float* Zs_out,
	float* VXs_out, float* VYs_out, float* VZs_out,
	float* AX_gravs_out, float* AY_gravs_out, float* AZ_gravs_out,
	float* AX_presss_out, float* AY_presss_out, float* AZ_presss_out,
	float* AX_viscs_out, float* AY_viscs_out, float* AZ_viscs_out,
	float* Masss_out, float* Densitys_out, float* Pressures_out,
	uint8_t* Types_out, float* BoundaryPsis_out,
	bool use_Z = false
);

void SampleStride(size_t grid, size_t block, size_t sample_size, size_t stride,
	const float* Xs, const float* Ys, const float* Zs,
	const float* VXs, const float* VYs, const float* VZs,
	const float* AX_gravs, const float* AY_gravs, const float* AZ_gravs,
	const float* AX_presss, const float* AY_presss, const float* AZ_presss,
	const float* AX_viscs, const float* AY_viscs, const float* AZ_viscs,
	const float* Masss, const float* Densitys, const float* Pressures,
	const uint8_t* Types, const float* BoundaryPsis,
	float* Xs_out, float* Ys_out, float* Zs_out,
	float* VXs_out, float* VYs_out, float* VZs_out,
	float* AX_gravs_out, float* AY_gravs_out, float* AZ_gravs_out,
	float* AX_presss_out, float* AY_presss_out, float* AZ_presss_out,
	float* AX_viscs_out, float* AY_viscs_out, float* AZ_viscs_out,
	float* Masss_out, float* Densitys_out, float* Pressures_out,
	uint8_t* Types_out, float* BoundaryPsis_out,
	bool use_Z = false);

#endif