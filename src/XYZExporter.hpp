#pragma once
#include "Observer.hpp"
#include <fstream>
#include <format>

class XYZExporter : public Observer<3> {
public:
	~XYZExporter() override = default;

	/// Executes on the simulation thread
	void OnEndFrame() override {
		if (m_Sim == nullptr)
			return;

		std::ofstream out(std::format("output-{}", m_Sim->GetFrame()));
		for (auto& p : m_Sim->GetParticles())
			out << p.Position.x << " " << p.Position.y << " " << p.Position.z << "\n";
	}

};
