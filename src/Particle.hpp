#pragma once
#include "Utility.hpp"

struct Particle
{
	using coord_t = coord<float, 2>;
	coord_t pos;
	float density;
};