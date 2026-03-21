#pragma once
#include "Base/Neighbors/NeighborFinder.hpp"


namespace mpi
{

template <size_t D, ParticleSet<D> Particles>
class SPHSimulation;

template <size_t D, ParticleSet<D> Particles>
class NeighborFinder : public base::NeighborFinder<D>
{
public:
	using idx_t = Particle<D>::idx_t;

	virtual ~NeighborFinder() override = default;

	virtual void InitializeFrame(SPHSimulation<D, Particles>* sim) {}

protected:

private:

};

}
