#include "SimInitializer.hpp"


template <size_t D>
class BoxInitializer;

template <>
class BoxInitializer<2> : public SimInitializer<2>
{
public:
	void Init(std::vector<Particle<2>>& out, float h) const override;
};

inline void BoxInitializer<2>::Init(std::vector<Particle<2>>& out, float h) const
{
	using vec_t = Particle<2>::vec_t;

	// Setup for "Dam breaking"
	int maxx = 100, maxy = 100;
	vec_t vertical_direction{ 0, -1 };
	vec_t zero_direction{ 0, 0 };
	for (float x = maxx/ 3; x < 2 * maxx / 3; x += 0.5)
	{
		for (float y = 25.3; y < 78; y += 0.5)
		{
			Particle<2> p;
			p.Type = FLUID;
			p.Position.y = float(y) / maxy;
			p.Position.x = float(x) / maxx;
			p.Velocity = zero_direction;
			// Arbitrary Mass to match with other constants
			p.Mass = 0.025f;
			p.A_grav = G_CONSTANT * vertical_direction;
			out.push_back(p);
		}
	}

	// Create multiple walls of the box to avoid leaks
	// TODO: Add walls options to properties
	for (int i = 0; i < 2; i++)
	{
		BuildWallAlongX(out, 0.25 - h  * i, 0, 1.0, 0.0025);
		BuildWallAlongX(out, 0.8 + h   * i, 0, 1.0, 0.0025);
		BuildWallAlongY(out, 0 + h   * i, 0.25, 0.8, 0.0025);
	}
		BuildWallAlongY(out, 1.0, 0.25, 0.8, 0.0025);
		BuildBox(out, coord<float, 2>{0.1, 0.25}, 0.1, 0.1, 0.0025);
		BuildBox(out, coord<float, 2>{0.1, 0.4}, 0.1, 0.1, 0.0025);
}
