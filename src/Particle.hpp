#pragma once
#include "Utility.hpp"

enum ParticleType : u8
{
	FLUID,
	SOLID
};

/*
TODO: 
AoS   : [{x, y, z, vx, vy, vz, ...}]             Array of particles
Hybrid: {[{x, y, z}], [{vx, vy, vz}], ...}       One array per coord
SoA   : {[x], [y], [z], [vx], [vy], [vz], ...}   One array per component
*/

template <size_t D>
struct Particle
{
	using vec_t = coord<float, D>;
	using idx_t = u32; // alias to specify that the int represents the index of a particle

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
