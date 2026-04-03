#pragma once
#include "Utility.hpp"
#include "Object.hpp"

template <size_t D, ParticleSet<2> Particles>
class Source;

template <ParticleSet<2> Particles>
class Source<2, Particles> : public Object<2, Particles>
{
public:
	using idx_t = Object<2, Particles>::idx_t;
	Source(const coord<float, 2> A_, const coord<float, 2> B_, Particles& out_, const bool is_right_) : Object<2, Particles>(A_, B_, out_, is_right_) {}
	~Source() override = default;

	void OnFrameStart(const float time) override;
	void SetParams(const float v_, const int part_count_, const float delay_, const float mass_, const bool is_gravity_, const float max_time_);

	ObjectType GetType() const override { return ObjectType::SOURCE; }

private:
	float v;
	int part_count;
	float mass;
	bool is_gravity;
	float delay = std::numeric_limits<float>::max();
	float max_time = std::numeric_limits<float>::max();
	float last_activation = 0;
};

template <ParticleSet<2> Particles>
inline void Source<2, Particles>::OnFrameStart(const float time)
{
	if (time > max_time)
		return;
	
	if (time - last_activation < delay)
		return;
	last_activation = time;

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
		if (is_gravity)
			p.A_grav = G_CONSTANT * vertical;
		else 
			p.A_grav = {0, 0};
		this->out.PushBack(p);
	}
}

template <ParticleSet<2> Particles>
inline void Source<2, Particles>::SetParams(const float v_, const int part_count_, const float delay_, const float mass_, const bool is_gravity_, const float max_time_)
{
	v = v_;
	part_count = part_count_;
	delay = delay_;
	mass = mass_;
	is_gravity = is_gravity_;
	max_time = max_time_;
}

