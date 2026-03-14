#ifdef HAS_CUDA
#include "SpatialHashing.cuh"

namespace cudasph
{

	void mallocGridGPU(GridGPU& g, int maxParticles, int totalCells)
	{
		cudaMalloc(&g.cellHash, sizeof(unsigned int) * maxParticles);
		cudaMalloc(&g.particleIndex, sizeof(unsigned int) * maxParticles);

		cudaMalloc(&g.cellStart, sizeof(int) * totalCells);
		cudaMalloc(&g.cellEnd, sizeof(int) * totalCells);
	}
	void freeGridGPU(GridGPU& g)
	{
		cudaFree(g.cellHash);
		cudaFree(g.particleIndex);
		cudaFree(g.cellStart);
		cudaFree(g.cellEnd);
	}


	__device__ int3 calcCell(float3 p, Grid3 grid)
	{
		int3 cell;

		cell.x = floorf((p.x - grid.origin.x) / grid.cellSize);
		cell.y = floorf((p.y - grid.origin.y) / grid.cellSize);
		cell.z = floorf((p.z - grid.origin.z) / grid.cellSize);

		return cell;
	}
	__device__ int2 calcCell(float2 p, Grid2 grid)
	{
		int2 cell;

		cell.x = floorf((p.x - grid.origin.x) / grid.cellSize);
		cell.y = floorf((p.y - grid.origin.y) / grid.cellSize);

		return cell;
	}

	__device__ int flattenCell(int3 c, Grid3 grid)
	{
		return (c.z * grid.resolution.y + c.y) * grid.resolution.x + c.x;
	}
	__device__ int flattenCell(int2 c, Grid2 grid)
	{
		return c.y * grid.resolution.x + c.x;
	}

	__global__ void computeHashes(
		int N,
		float* xs, float* ys, float* zs,
		Grid3 grid,
		unsigned int* cellHash,
		unsigned int* particleIndex)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= N) return;

		float3 pos = make_float3(xs[i], ys[i], zs[i]);

		int3 cell = calcCell(pos, grid);
		int hash = flattenCell(cell, grid);

		cellHash[i] = hash;
		particleIndex[i] = i;
	}
	__global__ void computeHashes(
		int N,
		float* xs, float* ys,
		Grid2 grid,
		unsigned int* cellHash,
		unsigned int* particleIndex)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= N) return;

		float2 pos = make_float2(xs[i], ys[i]);

		int2 cell = calcCell(pos, grid);
		int hash = flattenCell(cell, grid);

		cellHash[i] = hash;
		particleIndex[i] = i;
	}

	__global__ void reorderParticles(
		int N,
		unsigned int* particleIndex,
		float* xs, float* ys, float* zs,
		float* xsSorted, float* ysSorted, float* zsSorted)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= N) return;

		int sortedIndex = particleIndex[i];

		xsSorted[i] = xs[sortedIndex];
		ysSorted[i] = ys[sortedIndex];
		zsSorted[i] = zs[sortedIndex];
	}
	__global__ void reorderParticles(
		int N,
		unsigned int* particleIndex,
		float* xs, float* ys,
		float* xsSorted, float* ysSorted)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= N) return;

		int sortedIndex = particleIndex[i];

		xsSorted[i] = xs[sortedIndex];
		ysSorted[i] = ys[sortedIndex];
	}

	__global__ void findCellStartEnd(
		int N,
		const unsigned int* cellHash,
		int* cellStart,
		int* cellEnd)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= N) return;

		unsigned int hash = cellHash[i];
		if (hash >= N) return;

		// First particle of this hash sets the start
		if (i == 0 || hash != cellHash[i - 1])
		{
			cellStart[hash] = i;
		}

		// Last particle of this hash sets the end
		if (i == N - 1 || hash != cellHash[i + 1])
		{
			cellEnd[hash] = i + 1;  // exclusive end
		}
	}

	void resetGridCells(GridGPU& grid, int totalCells)
	{
		cudaMemset(grid.cellStart, -1, sizeof(int) * totalCells);
		cudaMemset(grid.cellEnd, -1, sizeof(int) * totalCells);
	}


	void buildGrid(ParticlesCuda<3>& particles, Grid3& grid, GridGPU& gridGPU)
	{
		int N = particles.Count;
		using vec_t = Particle<3>::vec_t;

		int block = 256;
		int gridSize = (N + block - 1) / block;

		resetGridCells(gridGPU, grid.totalCells);

		computeHashes<<<gridSize, block>>>(
			N,
			particles.Xs,
			particles.Ys,
			particles.Zs,
			grid,
			gridGPU.cellHash,
			gridGPU.particleIndex
		);
		cudaError_t err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("computeHashes: %s\n", cudaGetErrorString(err));

		thrust::sort_by_key(
			thrust::device,
			gridGPU.cellHash,
			gridGPU.cellHash + N,
			gridGPU.particleIndex
		);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("sort: %s\n", cudaGetErrorString(err));

		float* XsSorted; cudaMalloc(&XsSorted, N * sizeof(float));
		float* YsSorted; cudaMalloc(&YsSorted, N * sizeof(float));
		float* ZsSorted; cudaMalloc(&ZsSorted, N * sizeof(float));

		reorderParticles<<<gridSize, block>>>(
			N,
			gridGPU.particleIndex,
			particles.Xs,
			particles.Ys,
			particles.Zs,
			XsSorted,
			YsSorted,
			ZsSorted
		);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("reorderParticles: %s\n", cudaGetErrorString(err));

		particles.Xs = XsSorted;
		particles.Ys = YsSorted;
		particles.Zs = ZsSorted;

		findCellStartEnd<<<gridSize, block>>>(
			N,
			gridGPU.cellHash,
			gridGPU.cellStart,
			gridGPU.cellEnd
		);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("findCellStartEnd: %s\n", cudaGetErrorString(err));
	}

	__global__ void print(int N, int* s)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= N) return;

		printf("%i", s[i]);
	}
	void buildGrid(ParticlesCuda<2>& particles, Grid2& grid, GridGPU& gridGPU)
	{
		int N = particles.Count;
		using vec_t = Particle<2>::vec_t;

		int block = 256;
		int gridSize = (N + block - 1) / block;

		resetGridCells(gridGPU, grid.totalCells);

		computeHashes<<<gridSize, block>>>(
			N,
			particles.Xs,
			particles.Ys,
			grid,
			gridGPU.cellHash,
			gridGPU.particleIndex
		);
		cudaError_t err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("computeHashes: %s\n", cudaGetErrorString(err));

		thrust::sort_by_key(
			thrust::device,
			gridGPU.cellHash,
			gridGPU.cellHash + N,
			gridGPU.particleIndex
		);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("sort %s\n", cudaGetErrorString(err));

		float* XsSorted; cudaMalloc(&XsSorted, N * sizeof(float));
		float* YsSorted; cudaMalloc(&YsSorted, N * sizeof(float));

		reorderParticles<<<gridSize, block>>>(
			N,
			gridGPU.particleIndex,
			particles.Xs,
			particles.Ys,
			XsSorted,
			YsSorted
		);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("reorderParticles: %s\n", cudaGetErrorString(err));

		particles.Xs = XsSorted;
		particles.Ys = YsSorted;

		findCellStartEnd<<<gridSize, block>>>(
			N,
			gridGPU.cellHash,
			gridGPU.cellStart,
			gridGPU.cellEnd
		);
		err = cudaDeviceSynchronize();
		if (err != cudaSuccess)
			printf("findCellStartEnd: %s\n", cudaGetErrorString(err));
	}

}
#endif