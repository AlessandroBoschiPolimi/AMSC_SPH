#include "SPHSimulation.hpp"



SPHSimulation::SPHSimulation()
{
	for (int x = 0; x < 100; x++)
	{
		for (int y = 0; y < 100; y++)
		{
			Particle p;
			p.pos.y = float(y) / 100;
			p.pos.x = float(x) / 100;
			p.density = p.pos.x / 2 + 0.1;
			m_Particles.push_back(p);
		}
	}
}


void SPHSimulation::Start()
{
	for (int i = 0; i < 10000; i++)
		Step();
}
void SPHSimulation::Step()
{
	NotifyStartFrame();

	coord<float, 2, true> coord1, coord2;
	auto res = coord1 + coord2;
	// artificial workload
	std::this_thread::sleep_for(100ms);
	for (auto& p : m_Particles)
	{
		p.density = std::min(to<f32>(1.0), p.density + m_Increment);
	}

	m_Time += m_TimeStep;
	m_Frame++;
	NotifyEndFrame();
}
