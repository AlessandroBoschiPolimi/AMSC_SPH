#pragma once
#include "Utility.hpp"

enum ParticleType : u8
{
	FLUID,
	SOLID
};

struct Particle
{
	using vec_t = coord<float, 2>;
	using idx_t = u32;
	vec_t Position;
	vec_t Velocity;
	vec_t A_grav;
	vec_t A_visc;
	vec_t A_press;
	float Mass;
	float Density;
	float Pressure;
	std::vector<idx_t> Neighbors;
	ParticleType Type;
	float BoundaryPsi;
};
