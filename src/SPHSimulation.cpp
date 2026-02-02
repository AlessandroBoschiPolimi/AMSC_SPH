#include "SPHSimulation.hpp"

static constexpr float G_CONSTANT = 9.81f;

/* References:
 * [1]https://cg.informatik.uni-freiburg.de/publications/2012_SIGGRAPH_rigidFluidCoupling.pdfa
 * [2]https://cg.informatik.uni-freiburg.de/course_notes/sim_10_sph.pdf
 */

SPHSimulation::SPHSimulation() : W_Ker(m_Params.SmoothingLength)
{
	//Setup for ''Dam breaking''
	int maxx = 100, maxy = 100;
	vec_t vertical_direction{0, -1};
	vec_t zero_direction{0, 0};
	for (float x = 2; x < maxx / 2; x += 0.5)
	{
		for (float y = 28; y < 78 ; y += 0.5)
		{
			Particle p;
			p.Type = FLUID;
			p.Position.y = float(y) / maxy;
			p.Position.x = float(x) / maxx;
			p.Velocity = zero_direction;
			// Arbitrary Mass to match with other constants
			p.Mass = 0.5f;
			p.A_grav = G_CONSTANT * vertical_direction;
			m_Particles.push_back(p);
		}
	}
	
	// Create multiple walls of the box to avoid leaks
	// TODO: Add walls options to properties
	float pos1 = 25.0f;
	float pos2 = 80.0f;
	for (int i = 0; i < 5; i++)
	{
		BuildXWall(pos1 - m_Params.SmoothingLength / 4 * i, maxx, maxy, 0, maxx);
		BuildXWall(pos2 - m_Params.SmoothingLength / 4 * i, maxx, maxy, 0, maxx);
		BuildYWall(       m_Params.SmoothingLength / 4 * i, maxx, maxy, pos1, pos2);
		BuildYWall(maxx - m_Params.SmoothingLength / 4 * i, maxx, maxy, pos1, pos2);
	}
}
void SPHSimulation::BuildXWall(float y, float maxx, float maxy, float begin, float end)
{	
	/*
	 * Builds a wall of SOLID particles at position y
	 */
	for (float x = begin; x <= end; x+=0.25)
	{
		Particle p;
		p.Type = SOLID;
		p.Position.y = float(y) / maxy;
		p.Position.x = float(x) / maxx;
		p.Velocity = vec_t{0.0f, 0.0f};
		m_Particles.push_back(p);
	}
}

void SPHSimulation::BuildYWall(float x, float maxx, float maxy, float begin, float end)
{	
	/*
	 * Builds a wall of SOLID particles at position x
	 */
	for (float y = begin; y <= end; y+=0.25)
	{
		Particle p;
		p.Type = SOLID;
		p.Position.y = float(y) / maxy;
		p.Position.x = float(x) / maxx;
		p.Velocity = vec_t{0.0f, 0.0f};
		m_Particles.push_back(p);
	}
}

void SPHSimulation::Start()
{
	// TODO: Specify final timestep, and display in the ui the current progress (percentage and value)
	for (int i = 0; i < 400000; i++){
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
	// We have to clear Grid at each timestep to avoid cummulation of neighbors
	m_Grid.clear();
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

	/* for 3d
	neighborOffsets[27] = {
		{-1,-1,-1}, {-1,-1,0}, {-1,-1,1},
		{-1, 0,-1}, {-1, 0,0}, {-1, 0,1},
		{-1, 1,-1}, {-1, 1,0}, {-1, 1,1},

		{ 0,-1,-1}, { 0,-1,0}, { 0,-1,1},
		{ 0, 0,-1}, { 0, 0,0}, { 0, 0,1},
		{ 0, 1,-1}, { 0, 1,0}, { 0, 1,1},

		{ 1,-1,-1}, { 1,-1,0}, { 1,-1,1},
		{ 1, 0,-1}, { 1, 0,0}, { 1, 0,1},
		{ 1, 1,-1}, { 1, 1,0}, { 1, 1,1}
	};
	*/

	// TODO: maybe reserve out to the number of particles / number of cells.
	// Same as above with grid
	out.clear();
	float h = m_Params.SmoothingLength;
	float h2 = h * h;

	cell_pos_t base = GetCellPosition(m_Particles[i].Position, m_Params.SmoothingLength);
	for (const cell_pos_t& off : neighborOffsets)
	{
		cell_pos_t c = base + off;

		auto it = m_Grid.find(c);
		if (it == m_Grid.end()) continue;

		for (idx_t j : it->second) {
			if (j == i) continue;

			vec_t r = m_Particles[j].Position - m_Particles[i].Position;
			if (Dot(r, r) < h2)
				out.push_back(j);
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
	 * Additionally, we need to compute 'Mass' of boundary particles (ParticlePsi) 
	 * In the first step, we also need to initialize density
	 */
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == SOLID)
			ComputeBoundaryPsi(i);
	}
	if (step_num == 0)
	{
		for (int i = 0; i < m_Particles.size(); i++)
		{
			if (m_Particles[i].Type == FLUID)
				ComputeDensity(i);
		}
	}
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == FLUID) {
			ComputeAccelerationViscosity(i);
		}
	}
	for (int i = 0; i < m_Particles.size(); i++)
	{
		if (m_Particles[i].Type == FLUID) {
			UpdateVelocityInitial(i);
			UpdatePositionInitial(i);
		}
	}
}
void SPHSimulation::IterativePressure()
{
	/*
	 * Use simple scheme with splitting
	 * After computing initial forces and moving particles, compute dansity and pressure 
	   and move particles again.
	 */
	for (int i = 0; i < m_Particles.size(); i++) {
		if (m_Particles[i].Type == SOLID)
			continue;
		ComputeDensity(i);
		ComputePressure(i);
		ComputeAccelerationPressure(i);
	}
	
	if (m_Command.Type != Command::NONE)
		for (int i = 0; i < m_Particles.size(); i++)
			EvaluateCommand(i);

	for (int i = 0; i < m_Particles.size(); i++) {
		if (m_Particles[i].Type == SOLID)
			continue;
		UpdateVelocityIteration(i);
		UpdatePositionIteration(i);
	}
}

void SPHSimulation::ComputeBoundaryPsi(idx_t i)
{
	/*
	 * Computes 'mass' of the boundary particles used
	   to implement collisions
	 */
	float V = 0;
	for (auto &j : m_Particles[i].Neighbors)
	{
		if (m_Particles[j].Type == FLUID)
			V += W_Ker.GetValue(m_Particles[i], m_Particles[j]);
	}
	//Clamp the values in case the volume is too small
	m_Particles[i].BoundaryPsi =  (V > 1e2) ?  m_Params.RestDensity / V : 0;
}
void SPHSimulation::ComputeDensity(idx_t i)
{
	/*
	 * Simple function to compute density
	 * Takes initial value of density produced by itself
	 * For security, clamps the value in the end to avoid disappearing particles
	 */
	m_Particles[i].Density = m_Particles[i].Mass * W_Ker.GetValue(m_Particles[i], m_Particles[i]);
	m_Particles[i].Density = 0.0f;
	for (auto &j : m_Particles[i].Neighbors) {
		float W_ij = W_Ker.GetValue(m_Particles[i], m_Particles[j]);
		if (m_Particles[j].Type == FLUID)
		{
			m_Particles[i].Density += m_Particles[j].Mass * W_ij;
		}
		//Boundary handling
		else
		{
			m_Particles[i].Density += m_Particles[j].BoundaryPsi * W_ij;
		}
	}
	m_Particles[i].Density = std::max(m_Particles[i].Density, 0.1f * m_Params.RestDensity);
}

void SPHSimulation::ComputePressure(idx_t i)
{
	/*
	 * (Andrew) Tait equation
	 * Stiffness constant is user defined
	 */
	m_Particles[i].Pressure = m_Params.Stiffness *
				(std::pow(m_Particles[i].Density /
				m_Params.RestDensity, 7.0f) - 1);
}

void SPHSimulation::ComputeAccelerationViscosity(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [1]
	 * [2] has typos in that formula
	 */
	m_Particles[i].A_visc = vec_t{0, 0};
	for (auto &j : m_Particles[i].Neighbors)
	{
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		vec_t v_ij = m_Particles[i].Velocity - m_Particles[j].Velocity;
		vec_t x_ij = m_Particles[i].Position - m_Particles[j].Position;
		float prod = Dot(x_ij, v_ij);
		float numerator	= ((prod > 0) ? prod : 0) ;
		if (m_Particles[j].Type == FLUID)
		{
			numerator *= m_Particles[j].Mass;
			float denominator = (m_Particles[i].Density + m_Particles[j].Density) *(Dot(x_ij, x_ij) + 
						0.01 * m_Params.SmoothingLength * m_Params.SmoothingLength);
			m_Particles[i].A_visc +=  (2 * m_Params.Viscosity * (numerator / denominator)) * DW_ij;
		}
		else
		{
			float denominator = 2 * m_Particles[i].Density * (Dot(x_ij, x_ij) + 
						        0.01 * m_Params.SmoothingLength * m_Params.SmoothingLength);
			m_Particles[i].A_visc += (m_Params.ViscosityRigid * (numerator / denominator)) *
									  m_Particles[j].BoundaryPsi * DW_ij;
		}
	}
}

void SPHSimulation::ComputeAccelerationPressure(idx_t i)
{
	/*
	 * Computes acceleration dues to viscosity from [2]
	 * Boundary handling according to [2]
	 */
	m_Particles[i].A_press = vec_t{0, 0};
	for (auto &j : m_Particles[i].Neighbors)
	{
		float factor = (m_Particles[j].Type == SOLID) ? (m_Particles[j].BoundaryPsi) : m_Particles[j].Mass;
		idx_t z = (m_Particles[j].Type == SOLID) ? i : j;
		vec_t DW_ij = W_Ker.GetGradient(m_Particles[i], m_Particles[j]);
		m_Particles[i].A_press -= 
				(m_Particles[i].Pressure / (m_Particles[i].Density * m_Particles[i].Density) +
				m_Particles[z].Pressure / (m_Particles[z].Density * m_Particles[z].Density)) *
				factor *
				DW_ij;
	}
}

void SPHSimulation::UpdatePositionInitial(idx_t i)
{
	//Update position due to viscosity and gravity
	m_Particles[i].Position +=
				  m_Params.TimeStep *
				  m_Particles[i].Velocity;
}
void SPHSimulation::UpdatePositionIteration(idx_t i)
{
	//Update position due to pressure
	m_Particles[i].Position +=
				  m_Params.TimeStep *
				  m_Params.TimeStep *
				  m_Particles[i].A_press;
}
void SPHSimulation::UpdateVelocityInitial(idx_t i)
{
	//Update velocity due to viscosity and gravity
	m_Particles[i].Velocity += 
				  m_Params.TimeStep *
				  (m_Particles[i].A_grav +
				  m_Particles[i].A_visc);
}
void SPHSimulation::UpdateVelocityIteration(idx_t i)
{
	//Update velocity due to pressure
	m_Particles[i].Velocity +=
				  m_Params.TimeStep *
				  m_Particles[i].A_press;
}


void SPHSimulation::EvaluateCommand(idx_t i)
{
	float d = Norm(m_Particles[i].Position - m_Command.Position);
	if (d < m_Command.Radius)
	{
		float falloff = 1.0f - d / m_Command.Radius;
		m_Particles[i].Pressure += falloff * m_Command.Strength;
	}
}


SPHSimulation::cell_pos_t SPHSimulation::GetCellPosition(const vec_t& p, const float h)
{
	return {
		to<int>(std::floor(p.x / h)),
		to<int>(std::floor(p.y / h))
	};
}