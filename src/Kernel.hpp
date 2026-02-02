#pragma once
#include "Particle.hpp"
#include <cmath>
#include <numbers>


template <size_t D>
class Kernel {
public:
	using vec_t = Particle<D>::vec_t;

	Kernel(float h_);
	float GetValue(Particle<D>& xi, Particle<D>& xj) const;
	vec_t GetGradient(Particle<D>& xi, Particle<D>& xj) const;

private:
	// particle spacing
	float h;
	// Normalization constant to have integral over R^2 = 1
	float alpha;
};





template <size_t D>
Kernel<D>::Kernel(float h_) : h(h_)
{
	alpha = 5.0f / (14.0f * std::numbers::pi * h * h);
}

template <size_t D>
float Kernel<D>::GetValue(Particle<D>& xi, Particle<D>& xj) const
{
	/*
	Definition of the kernel function according
	to the presentation slide 52 (cubic slide)
	*/
	vec_t r = xi.Position - xj.Position;
	float q = std::sqrt(Dot(r, r)) / h;
	float out = 0.0f;
	if (q < 1.0f)
	{
		//faster then std::pow()
		float x1 = 2.0f - q;
		float x2 = 1.0f - q;
		out = x1 * x1 * x1 - 4.0f * x2 * x2 * x2;
	}
	else if (q >= 1.0f && q < 2.0f)
	{
		float x1 = 2.0f - q;
		out = x1 * x1 * x1;
	}
	return alpha * out;
}

template <size_t D>
Kernel<D>::vec_t Kernel<D>::GetGradient(Particle<D>& xi, Particle<D>& xj) const
{
	/*
	Gradient of the kernel function according
	to the presentation slide 53 (cubic slide)
	*/
	vec_t r = xi.Position - xj.Position;
	float q = std::sqrt(Dot(r, r)) / h;
	float dW_dq = 0.0f;
	if (q < 1.0f)
	{
		float x1 = 2.0f - q;
		float x2 = 1.0f - q;
		dW_dq = -3.0f * x1 * x1 + 12.0f * x2 * x2;
	}
	else if (q >= 1.0f && q < 2.0f)
	{
		float x1 = 2.0f - q;
		dW_dq = -3.0f * x1 * x1;
	}
	float prefac = alpha * dW_dq / (std::sqrt(Dot(r, r)) * h);
	vec_t grad_W = prefac * r;
	return grad_W;
}
