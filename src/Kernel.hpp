#pragma once
#include "Particle.hpp"

class Kernel{
	public:
		using vec_t = Particle::vec_t;
		Kernel(float h_);
		float GetValue(Particle &xi, Particle &xj) const;
		vec_t GetGradient(Particle &xi, Particle &xj) const;
	private:
		//particle spacing
		float h;
		//Normalization constant to have integral over R^2 = 1
		float alpha;
};
