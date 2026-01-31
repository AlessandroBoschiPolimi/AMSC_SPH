#include "Kernel.hpp"
#include <cmath>
Kernel::Kernel(float h_):
	h(h_)
{
	alpha = 5.0f / (14.0f * M_PI * std::pow(h, 2));
}

float Kernel::GetValue(Particle &xi, Particle &xj) const
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

coord<float, 2> Kernel::GetGradient(Particle &xi, Particle &xj) const
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
	float prefac =alpha *  dW_dq / (std::sqrt(Dot(r, r)) * h);
	vec_t grad_W = prefac * r;
	return grad_W;
}

