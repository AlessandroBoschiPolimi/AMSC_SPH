#pragma once
#include "Utility.hpp"

namespace base
{

struct SPHParams
{
	float RestDensity = 1000.0f;
	float Stiffness = 1e2f;
	float Viscosity = 1e-6f;
	float ViscosityRigid = 5e-3f;
	float TimeStep = 0.00017f;
	float SmoothingLength = 0.0035f;
	float PressureTol = 1e-2f;
	float FinalTime = 100.0;
};

struct SPHProfiling
{
	stdc::nanoseconds Neighbors = 0ns, Initialize = 0ns, IterativePressure = 0ns;
};

}
