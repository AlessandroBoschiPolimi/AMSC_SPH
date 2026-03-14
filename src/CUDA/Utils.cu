#ifdef HAS_CUDA
#include "Utils.cuh"

__global__ void SampleStrideKernel(
	size_t sample_size, size_t size, size_t stride,
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
	bool use_Z
) {
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= sample_size) return;
	size_t src_idx = idx * stride;
	if (src_idx >= size) return;

	Xs_out[idx] = Xs[src_idx];
	Ys_out[idx] = Ys[src_idx];
	if (use_Z) Zs_out[idx] = Zs[src_idx];

	VXs_out[idx] = VXs[src_idx];
	VYs_out[idx] = VYs[src_idx];
	if (use_Z) VZs_out[idx] = VZs[src_idx];

	AX_gravs_out[idx] = AX_gravs[src_idx];
	AY_gravs_out[idx] = AY_gravs[src_idx];
	if (use_Z) AZ_gravs_out[idx] = AZ_gravs[src_idx];

	AX_presss_out[idx] = AX_presss[src_idx];
	AY_presss_out[idx] = AY_presss[src_idx];
	if (use_Z) AZ_presss_out[idx] = AZ_presss[src_idx];

	AX_viscs_out[idx] = AX_viscs[src_idx];
	AY_viscs_out[idx] = AY_viscs[src_idx];
	if (use_Z) AZ_viscs_out[idx] = AZ_viscs[src_idx];

	Masss_out[idx] = Masss[src_idx];
	Densitys_out[idx] = Densitys[src_idx];
	Pressures_out[idx] = Pressures[src_idx];
	Types_out[idx] = Types[src_idx];
	BoundaryPsis_out[idx] = BoundaryPsis[src_idx];
}

void SampleStride(size_t grid, size_t block, size_t sample_size, size_t size, size_t stride,
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
	bool use_Z)
{
	SampleStrideKernel<<<grid, block>>>(
			sample_size, size, stride,
			Xs, Ys, Zs,
			VXs, VYs, VZs,
			AX_gravs, AY_gravs, AZ_gravs,
			AX_presss, AY_presss, AZ_presss,
			AX_viscs, AY_viscs, AZ_viscs,
			Masss, Densitys, Pressures,
			Types, BoundaryPsis,
			Xs_out, Ys_out, Zs_out,
			VXs_out, VYs_out, VZs_out,
			AX_gravs_out, AY_gravs_out, AZ_gravs_out,
			AX_presss_out, AY_presss_out, AZ_presss_out,
			AX_viscs_out, AY_viscs_out, AZ_viscs_out,
			Masss_out, Densitys_out, Pressures_out,
			Types_out, BoundaryPsis_out,
			use_Z
		);
}
#endif