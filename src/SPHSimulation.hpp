#pragma once
#include <vector>
#include <cmath>
#include "Observer.hpp"
#include "Kernel.hpp"

class SPHSimulation
{
public:
	using idx_t = u32;
	using cell_t = std::vector<idx_t>; // particle indices
	using cell_pos_t = coord<int, 2>;
	using grid_t = hmap<cell_pos_t, cell_t, CoordIntHash<2>>;
	using vec_t = Particle::vec_t;

	struct Params
	{
		float RestDensity = 1000.0f;
		float Stiffness = 0.001f;
		float Viscosity = 0.0003;
		float TimeStep = 0.001f;
		float SmoothingLength = 0.012f;
		float PressureTol = 1e-2f;
	};

public:
	SPHSimulation();
	~SPHSimulation()
	{
		for (auto* o : m_Observers)
			o->Attach(nullptr);
	}

	void Start();

	static cell_pos_t GetCellPosition(const vec_t& p, const float h);


private:
	void BuildYWall(float x, float maxx, float maxy, float begin, float end);
	void BuildXWall(float y, float maxx, float maxy, float begin, float end);

	void Step(int step_num);


	//Functions handling the three parts of the scheme for all the particles 
	void BuildGrid();
	void FindAllNeighbors();
	void Initialize(int step_num);
	void IterativePressure();
	//Helper functions for the subsequent particles
	void ComputeDensity(idx_t i);
	void ComputePressure(idx_t i);
	void ComputeAccelerationPressure(idx_t i);
	void ComputeAccelerationViscosity(idx_t i);
	void UpdatePositionInitial(idx_t i);
	void UpdatePositionIteration(idx_t i);
	void UpdateVelocityInitial(idx_t i);
	void UpdateVelocityIteration(idx_t i);
	/// Populates "out" with the neighbors of particle i-th
	void FindNeighbors(idx_t i, std::vector<idx_t>& out);


private:
	std::vector<Observer*> m_Observers;
	std::vector<Particle> m_Particles;
	size_t m_Frame = 0;
	float m_Time = 0.0f;

	Params m_Params;


	grid_t m_Grid;

	Kernel W_Ker;


public:
	void AddObserver(Observer* obs) {
		obs->Attach(this);
		m_Observers.push_back(obs);
	}

	void SetParams(const Params& params) { m_Params = params; }
	void SetRestDensity    (float val) { m_Params.RestDensity     = val; }
	void SetStiffness      (float val) { m_Params.Stiffness       = val; }
	void SetViscosity      (float val) { m_Params.Viscosity       = val; }
	void SetTimeStep       (float val) { m_Params.TimeStep        = val; }
	void SetSmoothingLength(float val) { m_Params.SmoothingLength = val; }

	Params GetParams() { return m_Params; }
	float GetRestDensity    () const { return m_Params.RestDensity    ; }
	float GetStiffness      () const { return m_Params.Stiffness      ; }
	float GetViscosity      () const { return m_Params.Viscosity      ; }
	float GetTimeStep       () const { return m_Params.TimeStep       ; }
	float GetSmoothingLength() const { return m_Params.SmoothingLength; }

	const std::vector<Particle>& GetParticles() const { return m_Particles; }
	const grid_t& GetGrid() const { return m_Grid; }

	float  GetTime()  const { return m_Time; }
	size_t GetFrame() const { return m_Frame; }


private:
	void NotifyStartFrame() {
		for (auto* o : m_Observers)
			o->OnStartFrame();
	}
	void NotifyEndFrame() {
		for (auto* o : m_Observers)
			o->OnEndFrame();
	}
};
