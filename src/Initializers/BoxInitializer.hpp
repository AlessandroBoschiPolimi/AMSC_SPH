#include "SimInitializer.hpp"


template <size_t D, ParticleSet<D> Particles>
class BoxInitializer;

template <ParticleSet<2> Particles>
class BoxInitializer<2, Particles> : public SimInitializer<2, Particles>
{
public:
	using Objects = SimInitializer<2, Particles>::Objects;

	void Init(Particles& out, float h, Objects& obj) const override;
};

template <ParticleSet<2> Particles>
inline void BoxInitializer<2, Particles>::Init(Particles& out, float h, Objects& obj) const
{
	using vec_t = Particle<2>::vec_t;

	// Setup for "Dam breaking"
	int maxx = 100, maxy = 100;
	vec_t vertical_direction{ 0, -1 };
	vec_t zero_direction{ 0, 0 };
	for (float x = 2; x < 84; x += 0.25)
	{
		for (float y = 5.3; y < 95; y += 0.25)
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
	// TODO: Add walls options to properties
	for (int i = 0; i < 2; i++)
	{
		this->BuildWallAlongX(out, 0.04 - h * i, 0, 1.0, 0.001);
		this->BuildWallAlongX(out, 0.98 + h * i, 0, 1.0, 0.001);
		this->BuildWallAlongY(out, 0 + h   * i, 0.04, 0.98, 0.001);
	}
	this->BuildWallAlongY(out, 1.0, 0.04, 0.98, 0.001);
	this->BuildBox(out, coord<float, 2>{0.85, 0.255}, 0.1, 0.1, 0.001);
	//this->BuildCircle(out, coord<float, 2>{0.92, 0.255}, 0.05, 0.001);

	//AddSink(obj, { 0.5, 0.2 }, { 0.7, 0.3 }, out, false);
	//AddSource(obj, { 0.4, 0.6 }, { 0.5, 0.5 }, out, true, 0.5, 10, 100, 0.025);
}
