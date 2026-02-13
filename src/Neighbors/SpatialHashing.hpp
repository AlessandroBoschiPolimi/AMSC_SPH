#pragma once
#include "NeighborFinder.hpp"


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
	const auto& particles = sim->GetParticles();
	float h = sim->GetSmoothingLength();
	
	#pragma omp master
	{
		m_Sim = sim;
		m_Grid.clear();
		m_Grid.reserve(particles.size());
	}
	#pragma omp barrier

	grid_t local_grid;
	local_grid.clear();
	local_grid.reserve(particles.size() / omp_get_num_threads());

	#pragma omp for
	for (int i = 0; i < particles.size(); i++) 
	{
		cell_pos_t c = GetCellPosition(particles[i].Position, h);
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

	const auto& particles = m_Sim->GetParticles();
	float h = m_Sim->GetSmoothingLength();
	float h2 = h * h;

	cell_pos_t base = GetCellPosition(particles[i].Position, h);
	for (const cell_pos_t& off : NeighborOffsets)
	{
		cell_pos_t c = base + off;

		auto it = m_Grid.find(c);
		if (it == m_Grid.end()) continue;

		for (idx_t j : it->second) {
			if (j == i) continue;

			vec_t r = particles[j].Position - particles[i].Position;
			if (Dot(r, r) < h2)
				out.push_back(j);
		}
	}
}


template <>
inline SpatialHashing<2>::cell_pos_t SpatialHashing<2>::GetCellPosition(const vec_t& p, const float h)
{
	return {
		to<int>(std::floor(p.x / h)),
		to<int>(std::floor(p.y / h))
	};
}
template <>
inline SpatialHashing<3>::cell_pos_t SpatialHashing<3>::GetCellPosition(const vec_t& p, const float h)
{
	return {
		to<int>(std::floor(p.x / h)),
		to<int>(std::floor(p.y / h)),
		to<int>(std::floor(p.z / h))
	};
}
