#pragma once
#include "Utility.hpp"

enum ParticleType : u8
{
	FLUID,
	SOLID
};

template <size_t D>
struct Particle
{
	using vec_t = coord<float, D>;
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
