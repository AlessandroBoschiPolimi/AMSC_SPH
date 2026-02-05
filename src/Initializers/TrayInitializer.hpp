#include "SimInitializer.hpp"


template <size_t D>
class TrayInitializer;

template <>
class TrayInitializer<3> : public SimInitializer<3>
{
public:
	void Init(std::vector<Particle<3>>& out, float h) const override;
};

inline void TrayInitializer<3>::Init(std::vector<Particle<3>>& out, float h) const
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
				for (float z = 5; z < maxz / 5; z += 0.5)
				{
					Particle<3> p;
					p.Type = FLUID;
					p.Position.x = float(x) / maxx;
					p.Position.y = float(y) / maxy;
					p.Position.z = float(z) / maxz;
					p.Velocity = zero_direction;
					// Arbitrary Mass to match with other constants
					p.Mass = 0.5f;
					p.A_grav = G_CONSTANT * vertical_direction;
					out.push_back(p);
				}
			}
		}
	}

	// Create multiple walls of the box to avoid leaks
	// TODO: Add walls options to properties
	{
		float maxx = 1.0 / 4, maxy = 1.0 / 4, maxz = 0.1;
		for (int i = 0; i < 5; i++)
		{
			BuildWallAlongXY(out, 0 + h / 400 * i, 0, maxx, 0, maxy, 0.0025);
			BuildWallAlongYZ(out, 0 + h / 400 * i, 0, maxx, 0, maxz, 0.0025);
			BuildWallAlongYZ(out, maxx - h / 400 * i, 0, maxx, 0, maxz, 0.0025);
			BuildWallAlongXZ(out, 0 + h / 400 * i, 0, maxy, 0, maxz, 0.0025);
			BuildWallAlongXZ(out, maxy - h / 400 * i, 0, maxy, 0, maxz, 0.0025);
		}
	}
}