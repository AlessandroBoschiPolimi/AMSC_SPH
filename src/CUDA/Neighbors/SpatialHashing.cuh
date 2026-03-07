#pragma once
//#include "NeighborFinder.cuh"
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <device_launch_parameters.h>
#include <cuda.h>

namespace cudasph
{
	
__host__ __device__ inline int3 operator+(const int3& a, const int3& b)
{
	return make_int3(a.x + b.x,
		a.y + b.y,
		a.z + b.z);
}


struct GridParams
{
    float cellSize;
    int3 gridSize;
    int numCells;
};
struct NeighborGrid
{
	int numParticles;

	GridParams grid;

	int* particleHash;
	int* particleIndex;

	int* cellStart;
	int* cellEnd;

	float* sortedXs, *sortedYs, *sortedZs;
	float* sortedVXs, * sortedVYs, * sortedVZs;
};

void initNeighborGrid(
	NeighborGrid& grid,
	int numParticles,
	int3 gridSize)
{
	grid.numParticles = numParticles;

	grid.grid.gridSize = gridSize;
	grid.grid.numCells = gridSize.x * gridSize.y * gridSize.z;

	cudaMalloc(&grid.particleHash, numParticles * sizeof(int));
	cudaMalloc(&grid.particleIndex, numParticles * sizeof(int));

	cudaMalloc(&grid.cellStart, grid.grid.numCells * sizeof(int));
	cudaMalloc(&grid.cellEnd, grid.grid.numCells * sizeof(int));

	cudaMalloc(&grid.sortedXs, numParticles * sizeof(float));
	cudaMalloc(&grid.sortedYs, numParticles * sizeof(float));
	cudaMalloc(&grid.sortedZs, numParticles * sizeof(float));
	cudaMalloc(&grid.sortedVXs, numParticles * sizeof(float));
	cudaMalloc(&grid.sortedVYs, numParticles * sizeof(float));
	cudaMalloc(&grid.sortedVZs, numParticles * sizeof(float));
}

__global__ void computeHashKernel(
	int n,
	const float* xs, const float* ys, const float* zs,
	int* particleHash,
	int* particleIndex,
	GridParams grid,
	int i);

__global__ void reorderParticlesKernel(
	int n,
	const int* particleIndex,
	const float* xs, const float* ys, const float* zs,
	const float* vxs, const float* vys, const float* vzs,
	float* sxs, float* sys, float* szs, // sorted
	float* svxs, float* svys, float* svzs,
	int i);

__global__ void findCellStartEndKernel(
	int n,
	const int* particleHash,
	int* cellStart,
	int* cellEnd,
	int i);

void buildNeighborGrid(
	NeighborGrid& grid,
	float3* pos,
	float3* vel);

}