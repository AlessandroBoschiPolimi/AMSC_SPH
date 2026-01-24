#pragma once
#include "Observer.hpp"
#include <fstream>
#include <format>

class XYZExporter : public Observer {
public:
	~XYZExporter() override = default;

	/// Executes on the simulation thread
	void OnEndFrame() override {
		// TODO: ASSERT m_Sim != nullptr

		std::ofstream out(std::format("output-{}", m_Sim->GetTime()));
		for (auto& p : *m_particles) {
			out << p.Position.x << " " << p.Position.y << " 0\n";
		}
	}


private:
	const std::vector<Particle>* m_particles = nullptr;
};
