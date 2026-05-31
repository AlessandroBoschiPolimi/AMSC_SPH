#pragma once
#include "SimInitializer.hpp"


template <size_t D, ParticleSet<D> Particles>
class PascalInitializer;

template <ParticleSet<2> Particles>
class PascalInitializer<2, Particles> : public SimInitializer<2, Particles>
{
public:
	using Objects = SimInitializer<2, Particles>::Objects;

	void Init(Particles& out, float h, Objects& obj) const override;

	std::pair<coord<float, 2>, coord<float, 2>> GetDomain() const override;
};

template <ParticleSet<2> Particles>
inline void PascalInitializer<2, Particles>::Init(Particles& out, float h, Objects& obj) const
{
	using vec_t = Particle<2>::vec_t;

	float density_wall = 0.002;
	this->AddSource(obj, { 0.11, 0.8501 }, { 0.18, 0.850 }, out, true, 1, 10, 0.0002f * 100, 0.01f, true, 0.0002f * 100000);
	
	this->BuildWallAlongY(out, 0.1 , 0.2, 0.85, density_wall);
	this->BuildWallAlongY(out, 0.2 , 0.25, 0.85, density_wall);
	this->BuildWallAlongY(out, 0.5 , 0.25, 0.85, density_wall);
	this->BuildWallAlongY(out, 0.8 , 0.2, 0.85, density_wall);
	
	this->BuildWallAlongX(out, 0.2 , 0.1, 0.8, density_wall);
	this->BuildWallAlongX(out, 0.25, 0.2, 0.5, density_wall);

}


template <ParticleSet<2> Particles>
inline std::pair<coord<float, 2>, coord<float, 2>> PascalInitializer<2, Particles>::GetDomain() const
{
	return std::pair<coord<float, 2>, coord<float, 2>>{{ 0.0f, 0.0f }, { 1.0f, 1.0f }};
}
