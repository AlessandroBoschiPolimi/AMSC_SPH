#pragma once
#include "Observer.hpp"
#include <fstream>
#include <format>


template <ParticleSet<3> Particles>
class XYZExporter : public Observer<3, Particles> {
public:
	~XYZExporter() override = default;

	/// Executes on the simulation thread
	void OnEndFrame() override;

private:
	stdc::time_point<stdclock> m_SimFrameEnd;
};



void writeXYZ(const int frame, const std::vector<Particle<3>>& particles);
void writeVTU(const int frame, const std::vector<Particle<3>>& particles);
void writeVTUBinary(const int frame, const std::vector<Particle<3>>& particles);



template <ParticleSet<3> Particles>
inline void XYZExporter<Particles>::OnEndFrame() {
	if (this->m_Sim == nullptr)
		return;

	std::cout << "Frame: " << this->m_Sim->GetFrame() << '\n';

	auto now = stdclock::now();
	std::cout << "Simulate ";
	print_time(now - m_SimFrameEnd);
	std::cout << ' ';

	if (this->m_Sim->GetFrame() % 10 == 0)
	{
		auto writeStart = stdclock::now();

		const auto& sim_particles = this->m_Sim->GetParticles();
		std::vector<Particle<3>> particles;
		particles.resize(sim_particles.Size());
		for (size_t i = 0; i < particles.size(); i++)
			particles[i] = sim_particles.Particle(i);

		writeVTUBinary(this->m_Sim->GetFrame() / 10, particles);
		
		auto writeEnd = stdclock::now();
		std::cout << "Write ";
		print_time(writeEnd - writeStart);
	}
	std::cout << '\n';

	m_SimFrameEnd = stdclock::now();
}
