#pragma once
#include "Particle.hpp"
#include <vector>
template <size_t D>
class Object;

template <>
class Object<2>
{
	public:
		using idx_t = u32;
		Object(const coord<float, 2> A_, const coord<float, 2> B_, std::vector<Particle<2>>& out_, const bool is_right_):
			out(out_),
			A(A_),
			B(B_),
			is_right(is_right_)
			{
				a = (A_.y - B_.y) / (A_.x - B_.x);
				b = A_.y - a*A_.x;
			}
		virtual void Activate(const idx_t i){};
		virtual void Activate(){};
	protected:
		const coord<float, 2> A;
		const coord<float, 2> B;
		float a;
		float b;
		const bool is_right;
		std::vector<Particle<2>>& out;
};
