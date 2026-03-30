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

template <size_t D, ParticleSet<D> Particles =ParticleAoS<D> >
class SpatialHashing : public NeighborFinder<D, Particles>
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

	void InitializeFrame(SPHSimulation<D, Particles>* sim) override;

	const grid_t& GetGrid() const { return m_Grid; }

	static cell_pos_t GetCellPosition(const vec_t& p, const float h);

private:
	grid_t m_Grid;
	SPHSimulation<D, Particles>* m_Sim = nullptr;

	static constexpr auto& NeighborOffsets = spatial_hashing_data_impl<D>::value;
};



template <size_t D, ParticleSet<D> Particles>
inline void SpatialHashing<D, Particles>::InitializeFrame(SPHSimulation<D, Particles>* sim)
{
	m_Sim = sim;

	const auto& particles_local = m_Sim->GetParticlesLocal();
	const auto& particles_ghost = m_Sim->GetParticlesGhost();
	int size = particles_local.Size() + particles_ghost.Size();
	float h = m_Sim->GetSmoothingLength();

	m_Grid.clear();
	m_Grid.reserve(size);

	for (size_t i = 0; i < size; ++i) {
		cell_pos_t c;
		if (i < particles_local.Size())
			c = GetCellPosition(particles_local.Position(i), h);
		else
			c = GetCellPosition(particles_ghost.Position(i - particles_local.Size()), h);
		m_Grid[c].push_back(i);
	}
}


template <size_t D, ParticleSet<D> Particles>
inline void SpatialHashing<D, Particles>::Find(idx_t i, std::vector<idx_t>& out)
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


template <size_t D, ParticleSet<D> Particles>
inline SpatialHashing<D, Particles>::cell_pos_t SpatialHashing<D, Particles>::GetCellPosition(const vec_t& p, const float h)
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
	else static_assert(false);
	return cell_pos_t{};
}


}
