#pragma once
#include <vector>
#include "Observer.hpp"


class SPHSimulation {
public:
	SPHSimulation();
	~SPHSimulation()
	{
		for (auto* o : m_Observers)
			o->Attach(nullptr);
	}

	void Start();

	void AddObserver(Observer* obs) {
		obs->Attach(this);
		m_Observers.push_back(obs);
	}

	void SetRestDensity(float val) { m_RestDensity = val; }
	void SetStiffness  (float val) { m_Stiffness   = val; }
	void SetViscosity  (float val) { m_Viscosity   = val; }
	void SetTimeStep   (float val) { m_TimeStep    = val; }
	void SetIncrement  (float val) { m_Increment   = val; }

	float GetRestDensity() const { return m_RestDensity; }
	float GetStiffness  () const { return m_Stiffness  ; }
	float GetViscosity  () const { return m_Viscosity  ; }
	float GetTimeStep   () const { return m_TimeStep   ; }
	float GetIncrement  () const { return m_Increment  ; }

	const std::vector<Particle>& GetParticles() const { return m_Particles; }

	float  GetTime () const { return m_Time ; }
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

	void Step();


private:
	std::vector<Observer*> m_Observers;
	std::vector<Particle> m_Particles;
	size_t m_Frame = 0;
	float m_Time = 0.0f;

	float m_RestDensity = 1000.0f;
	float m_Stiffness   = 2000.0f;
	float m_Viscosity   =    0.1f;
	float m_TimeStep    =    0.001f;
	float m_Increment   =    0;
};
