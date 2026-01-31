#include "SPHSimulation.hpp"

#define G_CONSTANT 9.81f


SPHSimulation::SPHSimulation():
	W_Ker(m_Params.SmoothingLength)
{
	int maxx = 100, maxy = 100;
	vec_t vertical_direction{0, -1};
	vec_t zero_direction{0, 0};
	for (int x = 2; x < maxx / 2; x++)
	{
		for (int y = 30; y < 60 ; y++)
		{
			Particle p;
			p.Type = FLUID;
			p.Position.y = float(y) / maxy;
			p.Position.x = float(x) / maxx;
			p.Velocity = zero_direction;
			p.Mass = 0.1f;
			p.A_grav = G_CONSTANT * vertical_direction;
			m_Particles.push_back(p);
		}
	}
	float pos1 = 20.0f;
	float pos2 = 80.0f;
//	for (int i = 0; i < 2; i++)
//	{
//		BuildXWall(pos1 -0.1 * i, maxx, maxy, 0, maxx);
//		BuildXWall(pos2 - 0.1 * i, maxx, maxy, 0, maxx);
//		BuildYWall(0.1 * i, maxx, maxy, pos1, pos2);
//		BuildYWall(maxx - 0.1 * i, maxx, maxy, pos1, pos2);
//	}
}
void SPHSimulation::BuildXWall(float y, float maxx, float maxy, float begin, float end)
{	
	/*
	 * Builds a wall of SOLID particles at position y
	 */
	for (int x = begin; x <= end; x++)
			{
					Particle p;
					p.Type = SOLID;
					p.Position.y = float(y) / maxy;
					p.Position.x = float(x) / maxx;
					p.Mass = 0.1f;
					m_Particles.push_back(p);
			}


}

void SPHSimulation::BuildYWall(float x, float maxx, float maxy, float begin, float end)
{	
	/*
	 * Builds a wall of SOLID particles at position x
	 */
	for (int y = begin; y <= end; y++)
			{
					Particle p;
					p.Type = SOLID;
					p.Position.y = float(y) / maxy;
					p.Position.x = float(x) / maxx;
					p.Mass = 0.1f;
					m_Particles.push_back(p);
			}


}

void SPHSimulation::Start()
{
	// TODO: Specify final timestep, and display in the ui the current progress (percentage and value)
	for (int i = 0; i < 10000; i++){
		Step(i);
	}
}
void SPHSimulation::Step(int step_num)
{
	NotifyStartFrame();

	// Split the particles in a grid of subdomains.
	BuildGrid();
	FindAllNeighbors();
	Initialize(step_num);
	IterativePressure();	

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

void SPHSimulation::FindAllNeighbors()
{
	for (int i = 0; i < m_Particles.size(); i++)
		FindNeighbors(i, m_Particles[i].Neighbors);
}

void SPHSimulation::Initialize(int step_num)
{
	/*
	 * Non-iterative part of the timestep
	 * In each timestep, set initial force due to Viscosity
	 * Apply this and gravity force to all the particles
	 * In the first step, we also need to initialize density
	 */
	if (step_num == 0)
	{
		for (int i = 0; i < m_Particles.size(); i++)
			ComputeDensity(i);
	}
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == FLUID){
			ComputeAccelerationViscosity(i);
			UpdateVelocityInitial(i);
			UpdatePositionInitial(i);
		}
	}
}
void SPHSimulation::IterativePressure()
{
	/*
	 * Use iterative SESPH scheme with splitting, Presentation p. 121
	 * Recompute pressure and density multiple times until 
	 * the || rho_old - rho_new ||_max < err (user defined error)
	 * Unce scheem converges, keep the updated position end velocity
	 */
	float error = 2 * m_Params.PressureTol;
	while (error > m_Params.PressureTol){
		error = 0;
		for (int i = 0; i < m_Particles.size(); i++){
			float old_dens = m_Particles[i].Density;
			ComputeDensity(i);
			float new_dens = m_Particles[i].Density;
			float curr_error = std::abs((new_dens - old_dens) / old_dens);
			if (curr_error > error)
				error = curr_error;
			ComputePressure(i);
			if (m_Particles[i].Type == FLUID){
				ComputeAccelerationPressure(i);
				UpdateVelocityIteration(i);
				UpdatePositionIteration(i);
			}
		}
	}
}

void SPHSimulation::ComputeDensity(idx_t i)
{
	m_Particles[i].Density = 0;
	for (auto &j: m_Particles[i].Neighbors){
		float W_ij = W_Ker.GetValue(m_Particles[i], m_Particles[j]);
		m_Particles[i].Density +=  m_Particles[j].Mass* W_ij;
	}
}

void SPHSimulation::ComputePressure(idx_t i)
{
	//Calculate pressure from state equation, Presentation p. 121
	//Stiffness constant is user defined
	m_Particles[i].Pressure = m_Params.Stiffness *
				(m_Particles[i].Density -
				m_Params.RestDensity);
}

void SPHSimulation::ComputeAccelerationViscosity(idx_t i)
{
	//Compute acceleration due to viscosity, Presentation p. 102
	m_Particles[i].A_visc = vec_t{0, 0};
	for (auto &j: m_Particles[i].Neighbors){
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		vec_t v_ij = m_Particles[i].Velocity - m_Particles[j].Velocity;
		vec_t x_ij = m_Particles[i].Position - m_Particles[j].Position;
		float numerator	= m_Particles[j].Mass * Dot(x_ij, DW_ij);
		float denominator = m_Particles[j].Density * (Dot(x_ij, x_ij) + 
					0.01 * m_Params.SmoothingLength * m_Params.SmoothingLength);
		m_Particles[i].A_visc = m_Particles[i].A_visc + 
					(2 * m_Params.Viscosity *
					(numerator / denominator)) * v_ij;

	}
}

void SPHSimulation::ComputeAccelerationPressure(idx_t i)
{
	//Compute acceleration due to pressure, Presentation p. 102
	m_Particles[i].A_press = vec_t{0, 0};
	for (auto &j: m_Particles[i].Neighbors)
	{
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		m_Particles[i].A_press = m_Particles[i].A_press -
				m_Particles[j].Mass * 
				(m_Particles[j].Pressure / (m_Particles[j].Density * m_Particles[j].Density) +
				m_Particles[i].Pressure / (m_Particles[i].Density * m_Particles[i].Density)) *
				DW_ij;
	}
}

void SPHSimulation::UpdatePositionInitial(idx_t i)
{
	//Update position due to viscosity and gravity
	m_Particles[i].Position = m_Particles[i].Position +
				  m_Params.TimeStep *
				  m_Particles[i].Velocity;
}
void SPHSimulation::UpdatePositionIteration(idx_t i)
{
	//Update position due to pressure
	m_Particles[i].Position = m_Particles[i].Position +
				  m_Params.TimeStep *
				  m_Params.TimeStep *
				  m_Particles[i].A_press;
}
void SPHSimulation::UpdateVelocityInitial(idx_t i)
{
	//Update velocity due to viscosity and gravity
	m_Particles[i].Velocity = m_Particles[i].Velocity +
				  m_Params.TimeStep *
				  (m_Particles[i].A_grav +
				  m_Particles[i].A_visc);
}
void SPHSimulation::UpdateVelocityIteration(idx_t i)
{
	//Update velocity due to pressure
	m_Particles[i].Velocity = m_Particles[i].Velocity +
				  m_Params.TimeStep *
				  m_Particles[i].A_press;
}


SPHSimulation::cell_pos_t SPHSimulation::GetCellPosition(const vec_t& p, const float h)
{
	return {
		to<int>(std::floor(p.x / h)),
		to<int>(std::floor(p.y / h))
	};
}
