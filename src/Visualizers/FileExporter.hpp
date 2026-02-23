#pragma once
#include "Observer.hpp"
#include <fstream>
#include <format>


template <size_t D, ParticleSet<D> Particles>
class FileExporter : public Observer<D, Particles> {
public:
	virtual ~FileExporter() override = default;

	/// Executes on the simulation thread
	void OnEndFrame() override;

private:
	stdc::time_point<stdclock> m_SimFrameEnd;
};



void writeXYZ(const int frame, const std::vector<Particle<3>>& particles);
void writeXYZ(const int frame, const std::vector<Particle<2>>& particles);
void writeVTU(const int frame, const std::vector<Particle<3>>& particles);
void writeVTU(const int frame, const std::vector<Particle<2>>& particles);
void writeVTUBinary(const int frame, const std::vector<Particle<3>>& particles);
void writeVTUBinary(const int frame, const std::vector<Particle<2>>& particles);



template <size_t D, ParticleSet<D> Particles>
inline void FileExporter<D, Particles>::OnEndFrame() {
	if (this->m_Sim == nullptr)
		return;

	u64 frame = this->m_Sim->GetFrame();
	std::cout << "Frame: " << frame << '\n';

	auto now = stdclock::now();
	std::cout << "Simulate ";
	print_time(now - m_SimFrameEnd);
	std::cout << ' ';

	if (frame % 10 == 0)
	{
		auto writeStart = stdclock::now();

		const auto& sim_particles = this->m_Sim->GetParticles();
		std::vector<Particle<D>> particles;
		this->m_Sim->GetParticles().GetParticles(particles);
		
		writeVTUBinary(frame / 10, particles);
		
		auto writeEnd = stdclock::now();
		std::cout << "Write ";
		print_time(writeEnd - writeStart);
	}
	std::cout << '\n';

	m_SimFrameEnd = stdclock::now();
}
