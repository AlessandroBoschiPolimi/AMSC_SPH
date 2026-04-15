#pragma once
#include "SimInitializer.hpp"


template <size_t D, ParticleSet<D> Particles>
class SmallBoxInitializer;

template <ParticleSet<2> Particles>
class SmallBoxInitializer<2, Particles> : public SimInitializer<2, Particles>
{
public:
	using Objects = SimInitializer<2, Particles>::Objects;

	void Init(Particles& out, float h, Objects& obj) const override;

	std::pair<coord<float, 2>, coord<float, 2>> GetDomain() const override;
};

template <ParticleSet<2> Particles>
inline void SmallBoxInitializer<2, Particles>::Init(Particles& out, float h, Objects& obj) const
{
	using vec_t = Particle<2>::vec_t;

	// Setup for "Dam breaking"
	int maxx = 100, maxy = 100;
	vec_t vertical_direction{ 0, -1 };
	vec_t zero_direction{ 0, 0 };
	for (float x = 20; x < 60; x += 0.25)
	{
		for (float y = 5.3; y < 50; y += 0.25)
		{
			Particle<2> p;
			p.Type = FLUID;
			p.Position.y = float(y) / maxy;
			p.Position.x = float(x) / maxx;
			p.Velocity = zero_direction;
			// Arbitrary Mass to match with other constants
			p.Mass = 0.01f;
			p.A_grav = G_CONSTANT * vertical_direction;
			out.PushBack(p);
		}
	}
	// Create multiple walls of the box to avoid leaks
	for (int i = 0; i < 2; i++)
	{
		this->BuildWallAlongX(out, 0.04 - h * i, 0.15, 0.65, 0.002);
		this->BuildWallAlongX(out, 0.98 + h * i, 0.15, 0.65, 0.002);
		this->BuildWallAlongY(out, 0.15 + h * i, 0.04, 0.98, 0.002);
	}
	this->BuildWallAlongY(out, 0.65, 0.04, 0.98, 0.002);

	this->AddSink(obj, { -10.0, 0.0 }, { 10.0, -10.0 }, out, false);
	this->AddSink(obj, { -10.0, 1.0 }, { 10.0, 11.0 }, out, true);
	this->AddSink(obj, { 0.0, -10.0 }, { -10.0, 10.0 }, out, false);
	this->AddSink(obj, { 1.0, -10.0 }, { 11.0, 10.0 }, out, false);
}


template <ParticleSet<2> Particles>
inline std::pair<coord<float, 2>, coord<float, 2>> SmallBoxInitializer<2, Particles>::GetDomain() const
{
	return std::pair<coord<float, 2>, coord<float, 2>>{{ 0.0f, 0.0f }, { 1.0f, 1.0f }};
}
