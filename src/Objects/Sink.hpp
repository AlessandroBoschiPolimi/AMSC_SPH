#pragma once
#include "Object.hpp"

template <size_t D>
class Sink;

template <>
class Sink<2> : public Object<2>
{
public:
	using idx_t = u32;
	using Object::Object;

	void Activate(const idx_t i) override;
};
inline void Sink<2>::Activate(const idx_t i)
{
	bool inside_x = out[i].Position.x < std::max(A.x, B.x) && out[i].Position.x > std::min(A.x, B.x);
	bool inside_y = out[i].Position.y < std::max(A.y, B.y) && out[i].Position.y > std::min(A.y, B.y);
	if (out[i].Type == FLUID)
	{
		if ( ((out[i].Position.y < a * out[i].Position.x + b) == is_right ) && inside_x && inside_y)
			out.erase(out.begin() + i);
	}
}

