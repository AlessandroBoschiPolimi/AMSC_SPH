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
		unsigned int* cellHash,
		int* cellStart,
		int* cellEnd)
	{
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		if (i >= N) return;

		unsigned int hash = cellHash[i];

		if (i == 0)
			cellStart[hash] = 0;
		else
		{
			unsigned int prevHash = cellHash[i - 1];

			if (hash != prevHash)
			{
				cellStart[hash] = i;
				cellEnd[prevHash] = i;
			}
		}

		if (i == N - 1)
			cellEnd[hash] = N;
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
		size_t sharedMem = sizeof(float3) * block;

		resetGridCells(gridGPU, grid.totalCells);

		computeHashes<<<gridSize, block, sharedMem>>>(
			N,
			particles.Xs,
			particles.Ys,
			particles.Zs,
			grid,
			gridGPU.cellHash,
			gridGPU.particleIndex
		);

		thrust::sort_by_key(
			thrust::device,
			gridGPU.cellHash,
			gridGPU.cellHash + N,
			gridGPU.particleIndex
		);

		reorderParticles << <gridSize, block, sharedMem >> > (
			N,
			gridGPU.particleIndex,
			particles.Xs,
			particles.Ys,
			particles.Zs,
			particles.Xs,
			particles.Ys,
			particles.Zs
		);

		findCellStartEnd << <gridSize, block, sharedMem >> > (
			N,
			gridGPU.cellHash,
			gridGPU.cellStart,
			gridGPU.cellEnd
		);
	}
	void buildGrid(ParticlesCuda<2>& particles, Grid2& grid, GridGPU& gridGPU)
	{
		int N = particles.Count;
		using vec_t = Particle<2>::vec_t;

		int block = 256;
		int gridSize = (N + block - 1) / block;
		size_t sharedMem = sizeof(float2) * block;

		resetGridCells(gridGPU, grid.totalCells);

		computeHashes<<<gridSize, block, sharedMem>>>(
			N,
			particles.Xs,
			particles.Ys,
			grid,
			gridGPU.cellHash,
			gridGPU.particleIndex
		);

		thrust::sort_by_key(
			thrust::device,
			gridGPU.cellHash,
			gridGPU.cellHash + N,
			gridGPU.particleIndex
		);

		reorderParticles<<<gridSize, block, sharedMem>>>(
			N,
			gridGPU.particleIndex,
			particles.Xs,
			particles.Ys,
			particles.Xs,
			particles.Ys
		);

		findCellStartEnd<<<gridSize, block, sharedMem>>>(
			N,
			gridGPU.cellHash,
			gridGPU.cellStart,
			gridGPU.cellEnd
		);
	}


//void initNeighborGrid(
//	NeighborGrid& grid,
//	int numParticles,
//	int3 gridSize)
//{
//	grid.numParticles = numParticles;
//
//	grid.grid.gridSize = gridSize;
//	grid.grid.numCells = gridSize.x * gridSize.y * gridSize.z;
//
//	cudaMalloc(&grid.particleHash, numParticles * sizeof(int));
//	cudaMalloc(&grid.particleIndex, numParticles * sizeof(int));
//
//	cudaMalloc(&grid.cellStart, grid.grid.numCells * sizeof(int));
//	cudaMalloc(&grid.cellEnd, grid.grid.numCells * sizeof(int));
//
//	cudaMalloc(&grid.sortedXs, numParticles * sizeof(float));
//	cudaMalloc(&grid.sortedYs, numParticles * sizeof(float));
//	cudaMalloc(&grid.sortedZs, numParticles * sizeof(float));
//	cudaMalloc(&grid.sortedVXs, numParticles * sizeof(float));
//	cudaMalloc(&grid.sortedVYs, numParticles * sizeof(float));
//	cudaMalloc(&grid.sortedVZs, numParticles * sizeof(float));
//}
//
//__device__ int3 calcGridPos(float3 p, float cellSize)
//{
//	return make_int3(
//		floorf(p.x / cellSize),
//		floorf(p.y / cellSize),
//		floorf(p.z / cellSize)
//	);
//}
//
//__device__ int calcGridHash(int3 gridPos, int3 gridSize)
//{
//	gridPos.x = (gridPos.x % gridSize.x + gridSize.x) % gridSize.x;
//	gridPos.y = (gridPos.y % gridSize.y + gridSize.y) % gridSize.y;
//	gridPos.z = (gridPos.z % gridSize.z + gridSize.z) % gridSize.z;
//
//	return gridPos.z * gridSize.y * gridSize.x +
//		gridPos.y * gridSize.x +
//		gridPos.x;
//}
//
//__global__ void computeHashKernel(
//	int n,
//	const float* xs, const float* ys, const float* zs,
//	int* particleHash,
//	int* particleIndex,
//	GridParams grid)
//{
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
//	if (i >= n) return;
//
//	float3 p = { xs[i], ys[i], zs[i] };
//
//	int3 gridPos = calcGridPos(p, grid.cellSize);
//	int hash = calcGridHash(gridPos, grid.gridSize);
//
//	particleHash[i] = hash;
//	particleIndex[i] = i;
//}
//
//__global__ void findCellStartEndKernel(
//	int n,
//	const int* particleHash,
//	int* cellStart,
//	int* cellEnd)
//{
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
//	if (i >= n) return;
//
//	int hash = particleHash[i];
//
//	if (i == 0)
//	{
//		cellStart[hash] = 0;
//	}
//	else
//	{
//		int prevHash = particleHash[i - 1];
//		if (hash != prevHash)
//		{
//			cellStart[hash] = i;
//			cellEnd[prevHash] = i;
//		}
//	}
//
//	if (i == n - 1)
//		cellEnd[hash] = n;
//}
//
//__global__ void reorderParticlesKernel(
//	int n,
//	const int* particleIndex,
//	const float* xs, const float* ys, const float* zs,
//	const float* vxs, const float* vys, const float* vzs,
//	float* sxs, float* sys, float* szs, // sorted
//	float* svxs, float* svys, float* svzs)
//{
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
//	if (i >= n) return;
//
//	int src = particleIndex[i];
//
//	sxs[i] = xs[src];
//	sys[i] = ys[src];
//	szs[i] = zs[src];
//	svxs[i] = vxs[src];
//	svys[i] = vys[src];
//	svzs[i] = vzs[src];
//}
//
//template <typename Callback>
//__device__ void forEachNeighborTiled(
//	int particleIndex,
//	const float* xs, const float* ys, const float* zs,
//	const int* cellStart,
//	const int* cellEnd,
//	GridParams grid,
//	float h,
//	float* sxs, float* sys, float* szs, // shared
//	int blockDimX, int threadIdxX,
//	Callback callback)
//{
//	float3 pos_i = { xs[particleIndex], ys[particleIndex], zs[particleIndex] };
//
//	int3 gridPos = calcGridPos(pos_i, grid.cellSize);
//
//	float h2 = h * h;
//
//	for (int z = -1; z <= 1; z++)
//	{
//		for (int y = -1; y <= 1; y++)
//		{
//			for (int x = -1; x <= 1; x++)
//			{
//				int3 neighborCell = gridPos + make_int3(x, y, z);
//
//				int hash = calcGridHash(neighborCell, grid.gridSize);
//
//				int start = cellStart[hash];
//				if (start == -1) continue;
//
//				int end = cellEnd[hash];
//
//				for (int tile = start; tile < end; tile += blockDimX)
//				{
//					int idx = tile + threadIdxX;
//
//					if (idx < end)
//					{
//						sxs[threadIdxX] = xs[idx];
//						sys[threadIdxX] = ys[idx];
//						szs[threadIdxX] = zs[idx];
//					}
//
//					__syncthreads();
//
//					int tileSize = min(blockDimX, end - tile);
//
//					for (int j = 0; j < tileSize; j++)
//					{
//						float3 pos_j = { sxs[j], sys[j], szs[j] };
//
//						float3 rij;
//						rij.x = pos_i.x - pos_j.x;
//						rij.y = pos_i.y - pos_j.y;
//						rij.z = pos_i.z - pos_j.z;
//
//						float r2 = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
//
//						if (r2 < h2 && r2 > 0)
//						{
//							callback(particleIndex, tile + j, rij, r2);
//						}
//					}
//
//					__syncthreads();
//				}
//			}
//		}
//	}
//}
//
//
//void buildNeighborGrid(
//	NeighborGrid& grid,
//	float* xs, float* ys, float* zs,
//	float* vxs, float* vys, float* vzs)
//{
//	int threads = 256;
//	int blocks = (grid.numParticles + threads - 1) / threads;
//
//	computeHashKernel<<<blocks, threads>>>(
//		grid.numParticles,
//		xs, ys, zs,
//		grid.particleHash,
//		grid.particleIndex,
//		grid.grid);
//
//	thrust::sort_by_key(
//		thrust::device,
//		grid.particleHash,
//		grid.particleHash + grid.numParticles,
//		grid.particleIndex);
//
//	reorderParticlesKernel<<<blocks, threads>>>(
//		grid.numParticles,
//		grid.particleIndex,
//		xs, ys, zs,
//		vxs, vys, vzs,
//		grid.sortedXs, grid.sortedYs, grid.sortedZs,
//		grid.sortedVXs, grid.sortedVYs, grid.sortedVZs);
//
//	cudaMemset(
//		grid.cellStart,
//		-1,
//		grid.grid.numCells * sizeof(int));
//
//	cudaMemset(
//		grid.cellEnd,
//		-1,
//		grid.grid.numCells * sizeof(int));
//
//	findCellStartEndKernel<<<blocks, threads>>>(
//		grid.numParticles,
//		grid.particleHash,
//		grid.cellStart,
//		grid.cellEnd);
//}

}
#endif