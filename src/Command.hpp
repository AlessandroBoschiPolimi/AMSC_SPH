#pragma once
#include "Particle.hpp"

template <size_t D>
struct Command {
	enum CommandType {
		NONE, PRESSURE
	};

	CommandType Type = NONE;
	Particle<D>::vec_t Position;
	float Radius;
	float Strength;
};