#pragma once
#include "NeighborFinder.hpp"

#include <cassert>
#include <algorithm>
#include <ranges>
#include <numeric>


namespace serial
{

template <size_t D>
struct morton_sorting_data_impl;
template <>
struct morton_sorting_data_impl<2> {
	inline static constexpr coord<int, 2> value[3 * 3] = {
		{ -1, -1 }, { -1,  0 }, { -1, +1 },
		{  0, -1 }, {  0,  0 }, {  0, +1 },
		{ +1, -1 }, { +1,  0 }, { +1, +1 }
	};
};
template <>
struct morton_sorting_data_impl<3> {
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


template <size_t D, ParticleSet<D> Particles>
class MortonSorting : public NeighborFinder<D, Particles>
{
public:
	using       idx_t = Particle<D>::idx_t;
	using       vec_t = Particle<D>::vec_t;
	using  cell_pos_t = coord<i32, D>;
	using cell_upos_t = coord<u32, D>;


public:
	MortonSorting(const vec_t& domain_min, const vec_t& domain_max) : m_DomainMin(domain_min), m_DomainMax(domain_max) {}
	virtual ~MortonSorting() override = default;

	void Find(idx_t i, std::vector<idx_t>& out) override;

	void InitializeFrame(SPHSimulation<D, Particles>* sim) override;

	static cell_pos_t GetCellPosition(const vec_t& p, const float h);
	static u64 ExpandBits(u32 v);
	static u64 MortonCode(const cell_upos_t& pos);


private:
	/// Requires m_MinGrid and m_MaxGrid to be up to date, so call only after InitializeFrame (for the current frame)
	bool OutsideDomain(const cell_pos_t& pos);


private:
	SPHSimulation<D, Particles>* m_Sim = nullptr;

	const vec_t m_DomainMin, m_DomainMax;
	cell_pos_t m_MinGrid, m_MaxGrid;

	std::vector<size_t> m_CellStart;
	std::vector<size_t> m_CellEnd;
	std::vector<uint64_t> m_UniqueCells;

	uint64_t m_MinCode;
	uint64_t m_MaxCode;
	std::vector<int> m_MortonLookup;

	static constexpr auto& NeighborOffsets = morton_sorting_data_impl<D>::value;
};


template <size_t D, ParticleSet<D> Particles>
inline void MortonSorting<D, Particles>::InitializeFrame(SPHSimulation<D, Particles>* sim)
{
	m_Sim = sim;

	Particles& particles = sim->GetParticles();
	float h = m_Sim->GetSmoothingLength();
	size_t size = particles.Size();

	// Get each particle Morton code
	std::vector<u64> morton;
	morton.resize(size);

	m_MinGrid = GetCellPosition(m_DomainMin, h);
	m_MaxGrid = GetCellPosition(m_DomainMax, h);

	for (size_t i = 0; i < size; ++i)
	{
		cell_pos_t cell = GetCellPosition(particles.Position(i), h);
		if (OutsideDomain(cell))
			continue;

		morton[i] = MortonCode(to<cell_upos_t>(cell - m_MinGrid));
	}

	// Sort by Morton code
	{
		std::vector<size_t> indices;
		indices.resize(size);
		std::iota(indices.begin(), indices.end(), 0);

		std::sort(indices.begin(), indices.end(),
			[&](size_t a, size_t b)
			{
				return morton[a] < morton[b];
			});

		{
			Particles tmp;
			tmp.Resize(size);
			for (size_t i = 0; i < indices.size(); ++i)
				tmp.SetParticle(i, particles, indices[i]);
			particles = std::move(tmp);
		}

		{
			std::vector<u64> tmp(size);
			for (size_t i = 0; i < indices.size(); ++i)
				tmp[i] = morton[indices[i]];
			morton = std::move(tmp);
		}
	}

	// Build cell ranges
	m_UniqueCells.clear();
	m_CellStart.clear();
	m_CellEnd.clear();
	for (size_t i = 0; i < size; )
	{
		u64 code = morton[i];
		size_t start = i;

		while (i < size && morton[i] == code)
			++i;

		size_t end = i;

		m_UniqueCells.push_back(code);
		m_CellStart.push_back(start);
		m_CellEnd.push_back(end);
	}

	m_MinCode = morton.front();
	m_MaxCode = morton.back();
	size_t tableSize = to<size_t>(m_MaxCode - m_MinCode + 1);

	// dense table: maps Morton code -> index in m_UniqueCells
	m_MortonLookup.clear();
	m_MortonLookup.resize(tableSize, -1);
	for (size_t idx = 0; idx < m_UniqueCells.size(); ++idx) {
		uint64_t code = m_UniqueCells[idx];
		m_MortonLookup[code - m_MinCode] = static_cast<int>(idx);
	}
}


template <size_t D, ParticleSet<D> Particles>
inline bool MortonSorting<D, Particles>::OutsideDomain(const cell_pos_t& pos)
{
	if constexpr (D == 2)
		return	pos.x < m_MinGrid.x || pos.x > m_MaxGrid.x ||
				pos.y < m_MinGrid.y || pos.y > m_MaxGrid.y;
	else if constexpr (D == 3)
		return	pos.x < m_MinGrid.x || pos.x > m_MaxGrid.x ||
				pos.y < m_MinGrid.y || pos.y > m_MaxGrid.y ||
				pos.z < m_MinGrid.z || pos.z > m_MaxGrid.z;
	else static_assert(D!=D);
	return false;
}


template <size_t D, ParticleSet<D> Particles>
inline void MortonSorting<D, Particles>::Find(idx_t i, std::vector<idx_t>& out)
{
	out.clear();

	const auto& particles = m_Sim->GetParticles();
	const float h = m_Sim->GetSmoothingLength();
	const float h2 = h * h;

	auto pi = particles.Position(i);
	cell_pos_t base = GetCellPosition(pi, h);

	for (const coord<int, D>& off : NeighborOffsets)
	{
		cell_pos_t neighbor = base + off;
		if (OutsideDomain(neighbor))
			continue;

		u64 neighborCode = MortonCode(to<cell_upos_t>(neighbor - m_MinGrid));

		// slower but safer when particles get out of bounds
		// with dense lookup map we just ignore oob particles (shouldn't we exterminate them anyway?)
		/*
		auto it = std::lower_bound(m_UniqueCells.begin(), m_UniqueCells.end(), neighborCode);
		if (it != m_UniqueCells.end() && *it == neighborCode)
		{
			size_t idx = std::distance(m_UniqueCells.begin(), it);

			for (size_t j = m_CellStart[idx]; j < m_CellEnd[idx]; ++j)
			{
				if (i == j) continue;
				
				auto r = particles[j].Position - particles[i].Position;
				if (Dot(r, r) < h2)
					out.push_back(j);
			}
		}
		*/

		if (neighborCode < m_MinCode || neighborCode > m_MaxCode)
			continue; // outside table

		int idx = m_MortonLookup[neighborCode - m_MinCode];
		if (idx == -1) continue; // empty cell

		// iterate particles in this neighbor cell
		for (size_t j = m_CellStart[idx]; j < m_CellEnd[idx]; ++j)
		{
			if (i == j) continue;
		
			auto r = particles.Position(j) - pi;
			if (Dot(r, r) < h2)
				out.push_back(j);
		}
	}
}

template <size_t D, ParticleSet<D> Particles>
inline MortonSorting<D, Particles>::cell_pos_t MortonSorting<D, Particles>::GetCellPosition(const vec_t& p, const float h)
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


template <size_t D, ParticleSet<D> Particles>
inline u64 MortonSorting<D, Particles>::ExpandBits(u32 v)
{
	if constexpr (D == 2)
	{
		v = (v | (v << 8)) & 0x00FF00FF;
		v = (v | (v << 4)) & 0x0F0F0F0F;
		v = (v | (v << 2)) & 0x33333333;
		v = (v | (v << 1)) & 0x55555555;
		return v;
	}
	else if constexpr (D == 3)
	{
		u64 x = v & 0x1fffff; // 21 bits
		x = (x | x << 32) & 0x1f00000000ffff;
		x = (x | x << 16) & 0x1f0000ff0000ff;
		x = (x | x << 8) & 0x100f00f00f00f00f;
		x = (x | x << 4) & 0x10c30c30c30c30c3;
		x = (x | x << 2) & 0x1249249249249249;
		return x;
	}
	else static_assert(D!=D);
	return v;
}
template <size_t D, ParticleSet<D> Particles>
inline u64 MortonSorting<D, Particles>::MortonCode(const cell_upos_t& pos)
{
	if constexpr (D == 2)
		return ExpandBits(pos.x) | (ExpandBits(pos.y) << 1);
	else if constexpr (D == 3)
		return ExpandBits(pos.x) | (ExpandBits(pos.y) << 1) | (ExpandBits(pos.z) << 2);
	return 0;
}


}