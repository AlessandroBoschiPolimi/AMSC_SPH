#pragma once
#include "Base/Neighbors/NeighborFinder.hpp"


namespace openmp
{

template <size_t D>
class SPHSimulation;

template <size_t D>
class NeighborFinder : public base::NeighborFinder<D>
{
public:
	using idx_t = Particle<D>::idx_t;

	virtual ~NeighborFinder() override = default;

	virtual void InitializeFrame(SPHSimulation<D>* sim) {}
	/// Populates "out" with the neighbors of particle i-th
	void Find(idx_t i, std::vector<idx_t>& out) = 0;

protected:

private:

};

}
