#pragma once
#ifdef HAS_CUDA
//#include "NeighborFinder.cuh"
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <device_launch_parameters.h>
#include <cuda.h>

#include "Particle.hpp"

namespace cudasph
{
	
__host__ __device__ inline int3 operator+(const int3& a, const int3& b)
{
	return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
__host__ __device__ inline int2 operator+(const int2& a, const int2& b)
{
	return make_int2(a.x + b.x, a.y + b.y);
}
__host__ __device__ inline Particle<3>::vec_t to_vec(const float3& a)
{
	return { a.x, a.y, a.z };
}
__host__ __device__ inline Particle<2>::vec_t to_vec(const float2& a)
{
	return { a.x, a.y };
}
template <typename T>
__device__ T dmin(const T& a, const T& b)
{
	return (a < b) ? a : b;
}

#ifndef __CUDACC__
// man i love ts
void __syncthreads();
#endif


struct Grid2
{
	float cellSize;
	float2 origin;

	int2 resolution;
	int totalCells;

	__host__ void init(float2 minBound, float2 maxBound, float smoothingLength)
	{
		cellSize = smoothingLength;
		origin = minBound;

		resolution.x = (int)ceil((maxBound.x - minBound.x) / cellSize);
		resolution.y = (int)ceil((maxBound.y - minBound.y) / cellSize);

		totalCells = resolution.x * resolution.y;
	}
};
struct Grid3
{
	float cellSize;
	float3 origin;

	int3 resolution;
	int totalCells;

	__host__ void init(float3 minBound, float3 maxBound, float smoothingLength)
	{
		cellSize = smoothingLength;
		origin = minBound;

		resolution.x = (int)ceil((maxBound.x - minBound.x) / cellSize);
		resolution.y = (int)ceil((maxBound.y - minBound.y) / cellSize);
		resolution.z = (int)ceil((maxBound.z - minBound.z) / cellSize);

		totalCells = resolution.x * resolution.y * resolution.z;
	}
};

struct GridGPU
{
	unsigned int* cellHash;
	unsigned int* particleIndex;

	int* cellStart;
	int* cellEnd;
};
void mallocGridGPU(GridGPU& g, int maxParticles, int totalCells);
void freeGridGPU(GridGPU& g);

__device__ int3 calcCell(float3 p, Grid3 grid);
__device__ int2 calcCell(float2 p, Grid2 grid);

__device__ int flattenCell(int3 c, Grid3 grid);
__device__ int flattenCell(int2 c, Grid2 grid);

__global__ void computeHashes(
	int N,
	float* xs, float* ys, float* zs,
	Grid3 grid,
	unsigned int* cellHash,
	unsigned int* particleIndex);
__global__ void computeHashes(
	int N,
	float* xs, float* ys,
	Grid2 grid,
	unsigned int* cellHash,
	unsigned int* particleIndex);

__global__ void reorderParticles(
	int N,
	unsigned int* particleIndex,
	float* xs, float* ys, float* zs,
	float* xsSorted, float* ysSorted, float* zsSorted);
__global__ void reorderParticles(
	int N,
	unsigned int* particleIndex,
	float* xs, float* ys,
	float* xsSorted, float* ysSorted);

__global__ void findCellStartEnd(
	int N,
	unsigned int* cellHash,
	int* cellStart,
	int* cellEnd);

void resetGridCells(GridGPU& grid, int totalCells);

void buildGrid(ParticlesCuda<3>& particles, Grid3& grid, GridGPU& gridGPU); 
void buildGrid(ParticlesCuda<2>& particles, Grid2& grid, GridGPU& gridGPU);


template <typename Callback>
__device__ void forEachNeighbor(
	int i,
	float3 pos_i,
	float h2,
	float* xs,
	float* ys,
	float* zs,
	int* cellStart,
	int* cellEnd,
	Grid3 grid,
	Callback callback)
{
	extern __shared__ float3 sharedPos3[];

	int tid = threadIdx.x;

	int3 cell = calcCell(pos_i, grid);

	for (int dz = -1; dz <= 1; dz++)
	for (int dy = -1; dy <= 1; dy++)
	for (int dx = -1; dx <= 1; dx++)
	{
		int3 neigh = make_int3(cell.x + dx, cell.y + dy, cell.z + dz);

		if (neigh.x < 0 || neigh.y < 0 || neigh.z < 0)
			continue;

		if (neigh.x >= grid.resolution.x ||
			neigh.y >= grid.resolution.y ||
			neigh.z >= grid.resolution.z)
			continue;

		int hash = flattenCell(neigh, grid);

		int start = cellStart[hash];
		int end   = cellEnd[hash];

		if (start == -1) continue;

		for (int j = start; j < end; j += blockDim.x)
		{
			int idx = j + tid;

			if (idx < end)
				sharedPos3[tid] = make_float3(xs[idx], ys[idx], zs[idx]);

			__syncthreads();

			int count = dmin((int)blockDim.x, end - j);

			#pragma unroll
			for (int k = 0; k < count; k++)
			{
				int neighborIndex = j + k;
				float3 pos_j = sharedPos3[k];

				float3 r;
				r.x = pos_i.x - pos_j.x;
				r.y = pos_i.y - pos_j.y;
				r.z = pos_i.z - pos_j.z;

				float r2 = r.x*r.x + r.y*r.y + r.z*r.z;

				if (r2 < h2)
				{
					callback(i, neighborIndex, to_vec(r), r2);
				}
			}

			__syncthreads();
		}
	}
}

template <typename Callback>
__device__ void forEachNeighbor(
	int i,
	float2 pos_i,
	float h2,
	float* xs,
	float* ys,
	int* cellStart,
	int* cellEnd,
	Grid2 grid,
	Callback callback)
{
	extern __shared__ float2 sharedPos2[];

	int tid = threadIdx.x;

	int2 cell = calcCell(pos_i, grid);

	for (int dy = -1; dy <= 1; dy++)
	for (int dx = -1; dx <= 1; dx++)
	{
		int2 neigh = make_int2(cell.x + dx, cell.y + dy);

		if (neigh.x < 0 || neigh.y < 0)
			continue;

		if (neigh.x >= grid.resolution.x ||
			neigh.y >= grid.resolution.y)
			continue;

		int hash = flattenCell(neigh, grid);

		int start = cellStart[hash];
		int end   = cellEnd[hash];

		if (start == -1) continue;

		for (int j = start; j < end; j += blockDim.x)
		{
			int idx = j + tid;

			if (idx < end)
				sharedPos2[tid] = make_float2(xs[idx], ys[idx]);

			__syncthreads();

			int count = dmin((int)blockDim.x, end - j);

			#pragma unroll
			for (int k = 0; k < count; k++)
			{
				int neighborIndex = j + k;
				float2 pos_j = sharedPos2[k];

				float2 r;
				r.x = pos_i.x - pos_j.x;
				r.y = pos_i.y - pos_j.y;

				float r2 = r.x*r.x + r.y*r.y;

				if (r2 < h2)
				{
					callback(i, neighborIndex, to_vec(r), r2);
				}
			}

			__syncthreads();
		}
	}
}

}
#endif