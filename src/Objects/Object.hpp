#pragma once
#include "Particle.hpp"
#include <vector>


enum class ObjectType
{
	UNKNOWN, SOURCE, SINK, PROBE
};

template <size_t D, ParticleSet<D> Particles>
class Object;

template <ParticleSet<2> Particles>
class Object<2, Particles>
{
public:
	using idx_t = Particle<2>::idx_t;

	Object(const coord<float, 2> A_, const coord<float, 2> B_, Particles& out_, const bool is_right_)
		: out(out_), A(A_), B(B_), is_right(is_right_)
	{
		a = (A_.y - B_.y) / (A_.x - B_.x);
		b = A_.y - a * A_.x;
	}
	virtual ~Object() = default;

	virtual void Activate(const idx_t i) {};
	virtual void OnFrameStart(const float time) {};

	virtual ObjectType GetType() const { return ObjectType::UNKNOWN; }
	std::pair<coord<float, 2>, coord<float, 2>> GetPosition() const { return { A, B }; }

protected:
	const coord<float, 2> A;
	const coord<float, 2> B;
	float a;
	float b;
	const bool is_right;
	
	Particles& out;
};


template <ParticleSet<3> Particles>
class Object<3, Particles>
{
public:
	using idx_t = Particle<3>::idx_t;

	Object() = default;
	virtual ~Object() = default;

	virtual void Activate(const idx_t i) {};
	virtual void OnFrameStart(const float time) {};

protected:
	Particles& out;
};
