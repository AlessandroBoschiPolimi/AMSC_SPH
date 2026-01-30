#include "SPHSimulation.hpp"

#define G_CONSTANT 9.81f


SPHSimulation::SPHSimulation():
	W_Ker(m_Params.SmoothingLength)
{
	int maxx = 100, maxy = 100;
	vec_t vertical_direction{0, 1};
	vec_t zero_direction{0, 0};
	for (int x = 0; x < maxx; x++)
	{
		for (int y = 0; y < maxy; y++)
		{
			Particle p;
			p.Position.y = float(y) / maxy;
			p.Position.x = float(x) / maxx;
			p.Density = p.Position.x / 2 + 0.1;
			p.Velocity = zero_direction;
			p.F_grav = G_CONSTANT * vertical_direction;
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
//	ComputeDensity();
//	ComputePressure();
//
//	ComputeForces();   // pressure + viscosity + gravity
//	Integrate();
//	HandleBoundaries();

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

void SPHSimulation::FindAllNeighbors(){
	for (int i = 0; i < m_Particles.size(); i++)
		FindNeighbors(i, m_Particles[i].Neighbors);
}

void SPHSimulation::ComputeDensity(idx_t i)
{
		m_Particles[i].Density = 0;
		for (auto &j: m_Particles[i].Neighbors){
			float W_ij = W_Ker.GetValue(m_Particles[i], m_Particles[j]);
			m_Particles[i].Density +=  m_Particles[j].Mass* W_ij;
		}
}

void SPHSimulation::ComputeForceViscosity(idx_t i)
{
	m_Particles[i].F_visc = vec_t{0, 0};
	for (auto &j: m_Particles[i].Neighbors){
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		vec_t v_ij = m_Particles[i].Velocity - m_Particles[j].Velocity;
		vec_t x_ij = m_Particles[i].Position - m_Particles[j].Position;
		float numerator	= m_Particles[j].Mass * Dot(x_ij, DW_ij);
		float denominator = m_Particles[j].Density * (Dot(x_ij, x_ij) + 
					0.01 * std::pow(m_Params.SmoothingLength, 2.0f));
		m_Particles[i].F_visc = m_Particles[i].F_visc + 
					2 * m_Params.Viscosity *
				       	m_Particles[i].Mass*	
					(numerator / denominator) * v_ij;
	}
}

void SPHSimulation::ComputeForcePressure(idx_t i)
{
	m_Particles[i].F_press = vec_t{0, 0};
	for (auto &j: m_Particles[i].Neighbors){
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		m_Particles[i].F_press = m_Particles[i].F_press -
				m_Particles[j].Mass * 
				(m_Particles[j].Pressure / std::pow(m_Particles[j].Density, 2.0f) +
				m_Particles[i].Pressure / std::pow(m_Particles[i].Density, 2.0f)) *
				DW_ij;
	}
	m_Particles[i].F_press = m_Particles[i].Mass /
				 m_Particles[i].Density *
				 m_Particles[i].F_press;
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
