#pragma once

namespace cudasph
{
	template <size_t D>
	class SPHSimulation;
}

template <size_t D>
class CudaObserver {
public:
	CudaObserver() = default;
	virtual ~CudaObserver() = default;

	/// Mainly used to read simulation parameters / timing
	/// Executes on the simulation thread
	virtual void OnStartFrame() {}
	/// Mainly used to alter simulation parameters / timing
	/// Executes on the simulation thread
	virtual void OnEndFrame() {}

	/// Executes on the simulation thread
	virtual void Attach(cudasph::SPHSimulation<D>* sim) { m_Sim = sim; }

protected:
	cudasph::SPHSimulation<D>* m_Sim = nullptr;
};
