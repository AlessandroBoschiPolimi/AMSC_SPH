#pragma once
#include "SimInitializer.hpp"


template <size_t D, ParticleSet<D> Particles>
class TrayInitializer;

template <ParticleSet<3> Particles>
class TrayInitializer<3, Particles> : public SimInitializer<3, Particles>
{
public:
	using Objects = SimInitializer<3, Particles>::Objects;

	void Init(Particles& out, float h, Objects& obj) const override;

	std::pair<coord<float, 3>, coord<float, 3>> GetDomain() const override;
};

template < ParticleSet<3> Particles>
inline void TrayInitializer<3, Particles>::Init(Particles& out, float h, Objects& obj) const
{
	using vec_t = Particle<3>::vec_t;

	// Setup for "Dam breaking"
	{
		int maxx = 100, maxy = 100, maxz = 100;
		vec_t vertical_direction{ 0, 0, -1 };
		vec_t zero_direction{ 0, 0, 0 };
		for (float x = 5; x < maxx / 5; x += 0.5)
		{
			for (float y = 5; y < maxy / 5; y += 0.5)
			{
				for (float z = 2; z < maxz / 5; z += 0.5)
				{
					Particle<3> p;
					p.Type = FLUID;
					p.Position.x = float(x) / maxx;
					p.Position.y = float(y) / maxy;
					p.Position.z = float(z) / maxz;
					p.Velocity = zero_direction;
					// Arbitrary Mass to match with other constants
					p.Mass = 0.0001f;
					p.A_grav = G_CONSTANT * vertical_direction;
					out.PushBack(p);
				}
			}
		}
	}

	// Create multiple walls of the box to avoid leaks
	{
		float maxx = 1.0 / 4, maxy = 1.0 / 4, maxz = 0.25;
		for (int i = 0; i < 1; i++)
		{
			this->BuildWallAlongXY(out, 0 + h  * i , 0, maxx, 0, maxy, 0.0025);
			this->BuildWallAlongYZ(out, 0 + h  * i, 0, maxx, 0, maxz, 0.0025);
			this->BuildWallAlongYZ(out, maxx - h * i, 0, maxx, 0, maxz, 0.0025);
			this->BuildWallAlongXZ(out, 0 + h * i, 0, maxy, 0, maxz, 0.0025);
			this->BuildWallAlongXZ(out, maxy - h  * i, 0, maxy, 0, maxz, 0.0025);
		}
	}
}

template< ParticleSet<3> Particles>
inline std::pair<coord<float, 3>, coord<float, 3>> TrayInitializer<3, Particles>::GetDomain() const
{
	return std::pair<coord<float, 3>, coord<float, 3>>{{ 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f }};
}
