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

	void OnFrameStart() override;
};
inline void Sink<2>::OnFrameStart()
{
	if (out.empty())
		return;

	std::stack<size_t> to_erase;

	for (size_t i = 0; i < out.size(); i++)
	{
		auto& p = out[i];

		bool inside_x = out[i].Position.x < std::max(A.x, B.x) && out[i].Position.x > std::min(A.x, B.x);
		bool inside_y = out[i].Position.y < std::max(A.y, B.y) && out[i].Position.y > std::min(A.y, B.y);

		if (out[i].Type == FLUID && ((out[i].Position.y < a * out[i].Position.x + b) == is_right) && inside_x && inside_y)
			to_erase.push(i);
	}

	size_t to_erase_count = to_erase.size();
	size_t read = out.size() - 1;
	while (!to_erase.empty())
	{
		size_t i = to_erase.top();
		to_erase.pop();

		out[i] = std::move(out[read--]);
	}

	out.resize(out.size() - to_erase_count);
}

