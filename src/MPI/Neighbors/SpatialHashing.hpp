#pragma once
#include "NeighborFinder.hpp"

namespace mpi
{

template <size_t D>
struct spatial_hashing_data_impl;
template <>
struct spatial_hashing_data_impl<2> {
	inline static constexpr coord<int, 2> value[3 * 3] = {
		{ -1, -1 }, { -1,  0 }, { -1, +1 },
		{  0, -1 }, {  0,  0 }, {  0, +1 },
		{ +1, -1 }, { +1,  0 }, { +1, +1 }
	};
};
template <>
struct spatial_hashing_data_impl<3> {
	inline static constexpr coord<int, 3> value[3 * 3 * 3] = {
		{-1,-1,-1}, {-1,-1,0}, {-1,-1,1},
		{-1, 0,-1}, {-1, 0,0}, {-1, 0,1},
		{-1, 1,-1}, {-1, 1,0}, {-1, 1,1},

		{ 0,-1,-1}, { 0,-1,0}, { 0,-1,1},
		{ 0, 0,-1}, { 0, 0,0}, { 0, 0,1},
		{ 0, 1,-1}, { 0, 1,0}, { 0, 1,1},

		{ 1,-1,-1}, { 1,-1,0}, { 1,-1,1},
		{ 1, 0,-1}, { 1, 0,0}, { 1, 0,1},
		{ 1, 1,-1}, { 1, 1,0}, { 1, 1,1}
	};
};

template <size_t D>
class SpatialHashing : public NeighborFinder<D>
{
public:
	using      idx_t = Particle<D>::idx_t;
	using      vec_t = Particle<D>::vec_t;
	using     cell_t = std::vector<idx_t>;
	using cell_pos_t = coord<int, D>;
	using     grid_t = hmap<cell_pos_t, cell_t, CoordIntHash<D>>;
	
	using Particles = ParticleAoS<D>;

public:
	virtual ~SpatialHashing() override = default;

	void Find(idx_t i, std::vector<idx_t>& out) override;

	void InitializeFrame(SPHSimulation<D>* sim) override;

	const grid_t& GetGrid() const { return m_Grid; }

	static cell_pos_t GetCellPosition(const vec_t& p, const float h);

private:
	grid_t m_Grid;
	SPHSimulation<D>* m_Sim = nullptr;

	static constexpr auto& NeighborOffsets = spatial_hashing_data_impl<D>::value;
};



template <size_t D>
inline void SpatialHashing<D>::InitializeFrame(SPHSimulation<D>* sim)
{
	const auto& particles_local = sim->GetParticlesLocal();
	const auto& particles_ghost = sim->GetParticlesGhost();
	float h = sim->GetSmoothingLength();
	size_t size = particles_local.Size() + particles_ghost.Size();

	#pragma omp master
	{
		m_Sim = sim;
		m_Grid.clear();
		m_Grid.reserve(size);
	}
	#pragma omp barrier

	grid_t local_grid;
	local_grid.clear();
	local_grid.reserve(size / omp_get_num_threads());

	#pragma omp for
	for (int i = 0; i < size; i++)
	{
		cell_pos_t c;
		if (i < particles_local.Size())
			c = GetCellPosition(particles_local.Position(i), h);
		else
			c = GetCellPosition(particles_ghost.Position(i - particles_local.Size()), h);
		local_grid[c].push_back(i);
	}

	#pragma omp critical
	{
		for (auto& [cell, vec] : local_grid)
			m_Grid[cell].insert(m_Grid[cell].end(), vec.begin(), vec.end());
	}
	#pragma omp barrier
}



template <size_t D>
inline void SpatialHashing<D>::Find(idx_t i, std::vector<idx_t>& out)
{
	// TODO: maybe reserve out to the number of particles / number of cells.
	out.clear();

	const auto& particles_local = m_Sim->GetParticlesLocal();
	const auto& particles_ghost = m_Sim->GetParticlesGhost();
	float h = m_Sim->GetSmoothingLength();
	float h2 = h * h;

	auto pi = particles_local.Position(i);
	cell_pos_t base = GetCellPosition(pi, h);
	for (const cell_pos_t& off : NeighborOffsets)
	{
		cell_pos_t c = base + off;

		auto it = m_Grid.find(c);
		if (it == m_Grid.end()) continue;

		for (idx_t j : it->second) {
			if (j == i) continue;
			vec_t r;
			if (j < particles_local.Size())
				r = particles_local.Position(j) - pi;
			else
				r = particles_ghost.Position(j - particles_local.Size()) - pi;
			if (Dot(r, r) < h2)
				out.push_back(j);
		}
	}
}


template <size_t D>
inline SpatialHashing<D>::cell_pos_t SpatialHashing<D>::GetCellPosition(const vec_t& p, const float h)
{
	if constexpr (D == 2)
	{
		return {
			to<i32>(std::floor(p.x / h)),
			to<i32>(std::floor(p.y / h))
		};
	}
	else if constexpr (D == 3)
	{
		return {
			to<i32>(std::floor(p.x / h)),
			to<i32>(std::floor(p.y / h)),
			to<i32>(std::floor(p.z / h))
		};
	}
	else static_assert(D!=D);
	return cell_pos_t{};
}


}
