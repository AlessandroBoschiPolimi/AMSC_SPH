#pragma once
#include "Utility.hpp"

struct Particle
{
	using vec_t = coord<float, 2>;
	using idx_t = u32;
	vec_t Position;
	vec_t Velocity;
	vec_t F_grav;
	vec_t F_visc;
	vec_t F_press;
	float Mass;
	float Density;
	float Pressure;
	std::vector<idx_t> Neighbors;
};
