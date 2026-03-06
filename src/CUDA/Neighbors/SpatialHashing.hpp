#pragma once
#include "NeighborFinder.hpp"
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

namespace cuda
{

struct SpatialHashGrid
{
	int numParticles;
	int numCells;

	float cellSize;

	u32* d_cellHash;
	u32* d_particleIndex;

	u32* d_cellStart;
	u32* d_cellEnd;

	float3* d_sortedPos;
};

__device__ inline int3 GetCellPosition(float3 p, float h)
{
	return make_int3(
		floorf(p.x / h),
		floorf(p.y / h),
		floorf(p.z / h)
	);
}

__device__ inline u32 HashCell(int3 c)
{
	const u32 p1 = 73856093;
	const u32 p2 = 19349663;
	const u32 p3 = 83492791;

	return (c.x * p1) ^ (c.y * p2) ^ (c.z * p3);
}

__global__
void ComputeHashesKernel(
	int N,
	float3* pos,
	float cellSize,
	u32* cellHash,
	u32* particleIndex)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= N) return;

	int3 cell = GetCellPosition(pos[i], cellSize);

	cellHash[i] = HashCell(cell);
	particleIndex[i] = i;
}

__global__
void BuildCellRangesKernel(
	int N,
	u32* cellHash,
	u32* cellStart,
	u32* cellEnd)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= N) return;

	u32 hash = cellHash[i];

	if (i == 0)
	{
		cellStart[hash] = 0;
	}
	else
	{
		u32 prevHash = cellHash[i - 1];

		if (prevHash != hash)
		{
			cellStart[hash] = i;
			cellEnd[prevHash] = i;
		}
	}

	if (i == N - 1)
	{
		cellEnd[hash] = N;
	}
}

__global__ void ReorderParticlesKernel(
	int N,
	float3* pos,
	float3* sortedPos,
	u32* particleIndex)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= N) return;

	u32 src = particleIndex[i];
	sortedPos[i] = pos[src];
}

__device__ const int3 NeighborOffsets[27] =
{
	{-1,-1,-1},{-1,-1,0},{-1,-1,1},
	{-1, 0,-1},{-1, 0,0},{-1, 0,1},
	{-1, 1,-1},{-1, 1,0},{-1, 1,1},

	{ 0,-1,-1},{ 0,-1,0},{ 0,-1,1},
	{ 0, 0,-1},{ 0, 0,0},{ 0, 0,1},
	{ 0, 1,-1},{ 0, 1,0},{ 0, 1,1},

	{ 1,-1,-1},{ 1,-1,0},{ 1,-1,1},
	{ 1, 0,-1},{ 1, 0,0},{ 1, 0,1},
	{ 1, 1,-1},{ 1, 1,0},{ 1, 1,1}
};

__global__ void NeighborSearchKernel(
	int N,
	float3* pos,
	float cellSize,
	float h2,

	u32* particleIndex,
	u32* cellStart,
	u32* cellEnd,

	int* neighborList,
	int* neighborCount,
	int maxNeighbors)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= N) return;

	float3 pi = pos[tid];
	int3 base = GetCellPosition(pi, cellSize);

	int count = 0;

	for (int n = 0; n < 27; n++)
	{
		int3 cell = make_int3(
			base.x + NeighborOffsets[n].x,
			base.y + NeighborOffsets[n].y,
			base.z + NeighborOffsets[n].z
		);

		u32 hash = HashCell(cell);

		int start = cellStart[hash];
		int end = cellEnd[hash];

		for (int j = start; j < end; j++)
		{
			int pj = particleIndex[j];

			if (pj == tid) continue;

			float3 r = pos[pj] - pi;

			float dist2 = r.x * r.x + r.y * r.y + r.z * r.z;

			if (dist2 < h2 && count < maxNeighbors)
			{
				neighborList[tid * maxNeighbors + count] = pj;
				count++;
			}
		}
	}

	neighborCount[tid] = count;
}

void BuildSpatialGrid(SpatialHashGrid& grid, float3* d_pos)
{
	int N = grid.numParticles;

	int threads = 256;
	int blocks = (N + threads - 1) / threads;

	ComputeHashesKernel<<<blocks, threads>>>(
			N,
			d_pos,
			grid.cellSize,
			grid.d_cellHash,
			grid.d_particleIndex
		);

	thrust::device_ptr<u32> hashPtr(grid.d_cellHash);
	thrust::device_ptr<u32> indexPtr(grid.d_particleIndex);

	thrust::sort_by_key(hashPtr, hashPtr + N, indexPtr);

	cudaMemset(grid.d_cellStart, 0xff, grid.numCells * sizeof(u32));
	cudaMemset(grid.d_cellEnd, 0, grid.numCells * sizeof(u32));

	BuildCellRangesKernel<<<blocks, threads>>>(
			N,
			grid.d_cellHash,
			grid.d_cellStart,
			grid.d_cellEnd
		);

	ReorderParticlesKernel<<<blocks, threads>>>(
			N,
			d_pos,
			grid.d_sortedPos,
			grid.d_particleIndex
		);
}

}