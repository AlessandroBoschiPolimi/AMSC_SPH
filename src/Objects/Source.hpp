#pragma once
#include "Utility.hpp"
#include "Object.hpp"

// TODO: make their behaviour depend on simulation time, not frame count

template <size_t D, ParticleSet<2> Particles>
class Source;

template <ParticleSet<2> Particles>
class Source<2, Particles> : public Object<2, Particles>
{
public:
	using idx_t = Object<2, Particles>::idx_t;

	~Source() override = default;

	void OnFrameStart() override;
	void SetParams(const float v_, const int part_count_, const int fp_count_, const float mass_);

private:
	float v;
	int part_count;
	int fp_count;
	float mass;
	int ii = 0;
};

template <ParticleSet<2> Particles>
inline void Source<2, Particles>::OnFrameStart()
{
	ii++;
	if (ii % fp_count != 0)
		return;
	ii = 0;
	const coord<float, 2> vertical{0, -1};
	const coord<float, 2> dlen = 1 / static_cast<float>(part_count) * (this->A - this->B);
	coord<float, 2> norm_v = v / std::sqrt(1 + this->a * this->a) * coord<float, 2>{-this->a, 1};
	if (!(this->is_right == (this->B + norm_v).y < this->a * (this->B + norm_v).x + this->b))
		norm_v = -1.0f * norm_v;
	for (int n = 0; n < part_count; n++)
	{
		Particle<2> p;
		const coord<float, 2> len = (this->B + (static_cast<float>(n) * dlen));
		p.Position.x = len.x;
		p.Position.y = len.y;
		p.Mass = mass;
		p.Velocity.x = norm_v.x;
		p.Velocity.y = norm_v.y;
		p.Density = p.Mass * Dot(len, len);
		p.Type = FLUID;
		p.A_grav = G_CONSTANT * vertical;
		this->out.PushBack(p);
	}
}

template <ParticleSet<2> Particles>
inline void Source<2, Particles>::SetParams(const float v_, const int part_count_, const int fp_count_, const float mass_)
{
	v = v_;
	part_count = part_count_;
	fp_count = fp_count_;
	mass = mass_;
}

