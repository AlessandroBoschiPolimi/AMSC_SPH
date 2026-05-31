#pragma once
#include "Observer.hpp"
#include <fstream>


enum ExportFormat
{
	VTU, VTU_01, XYZ, LANDMINE
};

template <size_t D, ParticleSet<D> Particles>
class FileExporter : public Observer<D, Particles> {
public:
	virtual ~FileExporter() override = default;

	/// Executes on the simulation thread
	void OnEndFrame() override;

	void SetFrequency(u32 freq) { m_Frequency = std::max(static_cast<u32>(1), freq); }
	void SetBaseName(const std::string& name) { m_BaseName = name; }
	void SetFormat(const ExportFormat& format) { m_Format = format; }

private:
	stdc::time_point<stdclock> m_SimFrameEnd;
	u32 m_Frequency = 100;
	std::string m_BaseName = "output";

	ExportFormat m_Format = VTU_01;
};



void writeXYZ(const std::string& basename, const int frame, const std::vector<Particle<3>>& particles);
void writeXYZ(const std::string& basename, const int frame, const std::vector<Particle<2>>& particles);
void writeVTU(const std::string& basename, const int frame, const std::vector<Particle<3>>& particles);
void writeVTU(const std::string& basename, const int frame, const std::vector<Particle<2>>& particles);
void writeVTUBinary(const std::string& basename, const int frame, const std::vector<Particle<3>>& particles);
void writeVTUBinary(const std::string& basename, const int frame, const std::vector<Particle<2>>& particles);
void writeForComparison(const std::string& basename, const int frame, const std::vector<Particle<3>>& particles);
void writeForComparison(const std::string& basename, const int frame, const std::vector<Particle<2>>& particles);



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

	if (frame % m_Frequency == 0)
	{
		auto writeStart = stdclock::now();

		const auto& sim_particles = this->m_Sim->GetParticles();
		std::vector<Particle<D>> particles;
		this->m_Sim->GetParticles().GetParticles(particles);
		
		switch (m_Format)
		{
		case VTU:
			writeVTU(m_BaseName, frame / m_Frequency, particles);
			break;
		case VTU_01:
			writeVTUBinary(m_BaseName, frame / m_Frequency, particles);
			break;
		case XYZ:
			writeXYZ(m_BaseName, frame / m_Frequency, particles);
			break;
		case LANDMINE:
			writeForComparison(m_BaseName, frame / m_Frequency, particles);
			break;
		}
		
		auto writeEnd = stdclock::now();
		std::cout << "Write ";
		print_time(writeEnd - writeStart);
	}
	std::cout << '\n';

	m_SimFrameEnd = stdclock::now();
}
