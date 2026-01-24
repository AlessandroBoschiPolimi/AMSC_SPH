#pragma once
#include "Utility.hpp"

struct Particle
{
	using vec_t = coord<float, 2>;
	vec_t Position;
	vec_t Velocity;
	vec_t Acceleration;
	float Mass;
	float Density;
	float Pressure;
};