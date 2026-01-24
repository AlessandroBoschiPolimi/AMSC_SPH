#include "SPHSimulation.hpp"



SPHSimulation::SPHSimulation()
{
	int maxx = 100, maxy = 100;
	for (int x = 0; x < maxx; x++)
	{
		for (int y = 0; y < maxy; y++)
		{
			Particle p;
			p.Position.y = float(y) / maxy;
			p.Position.x = float(x) / maxx;
			p.Density = p.Position.x / 2 + 0.1;
			m_Particles.push_back(p);
		}
	}
}


void SPHSimulation::Start()
{
	// TODO: Specify final timestep, and display in the ui the current progress (percentage and value)
	for (int i = 0; i < 10000; i++)
		Step();
}
void SPHSimulation::Step()
{
	NotifyStartFrame();

	// Split the particles in a grid of subdomains.
	BuildGrid();

	ComputeDensity();
	ComputePressure();

	ComputeForces();   // pressure + viscosity + gravity
	Integrate();
	HandleBoundaries();

	m_Time += m_Params.TimeStep;
	m_Frame++;
	NotifyEndFrame();
}


void SPHSimulation::BuildGrid()
{
	m_Grid.reserve(m_Particles.size());

	float h = m_Params.SmoothingLength;

	for (size_t i = 0; i < m_Particles.size(); ++i) {
		cell_pos_t c = GetCellPosition(m_Particles[i].Position, m_Params.SmoothingLength);
		m_Grid[c].push_back(i);
	}
}
void SPHSimulation::FindNeighbors(idx_t i, std::vector<idx_t>& out)
{
	static const cell_pos_t neighborOffsets[3 * 3] = {
		{ -1, -1 }, { -1,  0 }, { -1, +1 },
		{  0, -1 }, {  0,  0 }, {  0, +1 },
		{ +1, -1 }, { +1,  0 }, { +1, +1 }
	};

	// TODO: maybe reserve out to the number of particles / number of cells.

	float h = m_Params.SmoothingLength;
	float h2 = h * h;

	cell_pos_t base = GetCellPosition(m_Particles[i].Position, m_Params.SmoothingLength);

	for (const cell_pos_t& off : neighborOffsets) {
		cell_pos_t c = base + off;

		auto it = m_Grid.find(c);
		if (it == m_Grid.end()) continue;

		for (idx_t j : it->second) {
			if (j == i) continue;

			vec_t r = m_Particles[j].Position - m_Particles[i].Position;
			if (Dot(r, r) < h2) {
				out.push_back(j);
			}
		}
	}
}

void SPHSimulation::ComputeDensity()
{
}
void SPHSimulation::ComputePressure()
{
}

void SPHSimulation::ComputeForces()
{
}

void SPHSimulation::Integrate()
{
}

void SPHSimulation::HandleBoundaries()
{
}




SPHSimulation::cell_pos_t SPHSimulation::GetCellPosition(const vec_t& p, const float h)
{
	return {
		to<int>(std::floor(p.x / h)),
		to<int>(std::floor(p.y / h))
	};
}