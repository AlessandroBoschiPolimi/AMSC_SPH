#pragma once
#include "SimInitializer.hpp"


template <size_t D, ParticleSet<D> Particles>
class PipeInitializer;

template <ParticleSet<2> Particles>
class PipeInitializer<2, Particles> : public SimInitializer<2, Particles>
{
public:
	using Objects = SimInitializer<2, Particles>::Objects;

	void Init(Particles& out, float h, Objects& obj) const override;

	std::pair<coord<float, 2>, coord<float, 2>> GetDomain() const override;
};

template <ParticleSet<2> Particles>
inline void PipeInitializer<2, Particles>::Init(Particles& out, float h, Objects& obj) const
{
	using vec_t = Particle<2>::vec_t;

	float density_wall = 0.002;
	this->AddSink(obj, { 0.9, 0.5 }, { 0.95, 0.55 }, out, false);
	this->AddSource(obj, { 0.10, 0.10 }, { 0.1001, 0.15 }, out, true, 2, 10, 0.0002f * 30, 0.025, false);
	
	this->BuildWallAlongY(out, 0.085 , 0.1, 0.15, density_wall);
	this->BuildWallAlongX(out, 0.1 , 0.085, 0.8, density_wall);
	this->BuildWallAlongX(out, 0.15 , 0.085, 0.75, density_wall);

	this->BuildWallAlongY(out, 0.75 , 0.15, 0.25, density_wall);
	this->BuildWallAlongY(out, 0.8 , 0.1, 0.3, density_wall);

	this->BuildWallAlongX(out, 0.25 , 0.2, 0.75, density_wall);
	this->BuildWallAlongX(out, 0.3 , 0.25, 0.8, density_wall);

	this->BuildWallAlongY(out, 0.25 , 0.3, 0.45, density_wall);
	this->BuildWallAlongY(out, 0.2 , 0.25, 0.5, density_wall);
	
	this->BuildWallAlongX(out, 0.45 , 0.25, 0.35, density_wall);
	this->BuildWallAlongX(out, 0.5 , 0.2, 0.3, density_wall);
	
	this->BuildWallAlongY(out, 0.3 , 0.5, 0.9, density_wall);
	this->BuildWallAlongY(out, 0.35 , 0.45, 0.85, density_wall);

	this->BuildWallAlongX(out, 0.85 , 0.35, 0.6, density_wall);
	this->BuildWallAlongX(out, 0.9 , 0.3, 0.65, density_wall);
	
	this->BuildWallAlongY(out, 0.6 , 0.5, 0.85, density_wall);
	this->BuildWallAlongY(out, 0.65 , 0.55, 0.9, density_wall);
	
	this->BuildWallAlongX(out, 0.55 , 0.65, 0.95, density_wall);
	this->BuildWallAlongX(out, 0.5 , 0.6, 0.95, density_wall);

}


template <ParticleSet<2> Particles>
inline std::pair<coord<float, 2>, coord<float, 2>> PipeInitializer<2, Particles>::GetDomain() const
{
	return std::pair<coord<float, 2>, coord<float, 2>>{{ 0.0f, 0.0f }, { 1.0f, 1.0f }};
}
