#pragma once
#include "Particle.hpp"

/// Find hot neighbors in less than 100 meters from your location!

template <size_t D>
class NeighborFinder
{
public:
	using idx_t = Particle<D>::idx_t;

	virtual ~NeighborFinder() = default;

	virtual void InitializeFrame(SPHSimulation<D>* sim) {}
	/// Populates "out" with the neighbors of particle i-th
	virtual void Find(idx_t i, std::vector<idx_t>& out) = 0;

protected:

private:

};

