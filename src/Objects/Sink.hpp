#pragma once
#include "Object.hpp"

template <size_t D, ParticleSet<D> Particles>
class Sink;

template <ParticleSet<2> Particles>
class Sink<2, Particles> : public Object<2, Particles>
{
public:
	using idx_t = Object<2, Particles>::idx_t;
	
	Sink(const coord<float, 2> A_, const coord<float, 2> B_, Particles& out_, const bool is_right_) : Object<2, Particles>(A_, B_, out_, is_right_) {}
	~Sink() override = default;

	void OnFrameStart() override;
};

template <ParticleSet<2> Particles>
inline void Sink<2, Particles>::OnFrameStart()
{
	if (this->out.Empty())
		return;

	size_t write_end = this->out.Size();

	float mx = std::min(this->A.x, this->B.x), Mx = std::max(this->A.x, this->B.x);
	float my = std::min(this->A.y, this->B.y), My = std::max(this->A.y, this->B.y);

	for (size_t i = 0; i < write_end;)
	{
		auto x = this->out.PositionX(i);
		auto y = this->out.PositionY(i);

		bool inside_x = x < Mx && x > mx;
		bool inside_y = y < My && y > my;

		// TODO: render objects
		if (this->out.Type(i) == FLUID && ((y < this->a * x + this->b) == this->is_right) && inside_x && inside_y)
			this->out.SetParticle(i, std::move(this->out.GetParticle(--write_end))); // TODO: double check the move here
		else
			i++;
	}

	this->out.Resize(write_end);
}

