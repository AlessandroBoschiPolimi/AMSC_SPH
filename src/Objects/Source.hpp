#pragma once
#define G_CONSTANT 9.81f
#include "Object.hpp"

template<size_t D>
class Source;

template<>
class Source<2>: public Object<2>
{
	public:
		using idx_t = u32;
		void Activate() override;
		void SetParams(const float v_, const int part_count_, const int fp_count_, const float mass_);
		using Object::Object;
	private:
		float v;
		int part_count;
		int fp_count;
		float mass;
		int ii = 0;

};

inline void Source<2>::Activate()
{
	ii++;
	if (ii % fp_count != 0)
		return;
	ii = 0;
	const coord<float, 2> vertical{0, -1};
	const coord<float, 2> dlen = 1 / static_cast<float>(part_count) *
					(A - B);
	coord<float, 2> norm_v = v / std::sqrt(1 + a * a) * coord<float, 2>{-a, 1};
	if (!(is_right == (B + norm_v).y < a * (B + norm_v).x + b))
		norm_v = -1.0f * norm_v;
	for (int n = 0; n < part_count; n++)
	{
		Particle<2> p;
		const coord<float, 2> len = (B + (static_cast<float>(n) * dlen));
		p.Position.x = len.x;
		p.Position.y = len.y;
		p.Mass = mass;
		p.Velocity.x = norm_v.x;
		p.Velocity.y = norm_v.y;
		p.Density = p.Mass * Dot(len, len);
		p.Type = FLUID;
		p.A_grav = G_CONSTANT * vertical;
		out.push_back(p);
	}
}

inline void Source<2>::SetParams(const float v_, const int part_count_, const int fp_count_, const float mass_)
{
	v = v_;
	part_count = part_count_;
	fp_count = fp_count_;
	mass = mass_;
}

