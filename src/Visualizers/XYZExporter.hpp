#pragma once
#include "Observer.hpp"
#include <fstream>
#include <format>



class XYZExporter : public Observer<3> {
public:
	~XYZExporter() override = default;

	/// Executes on the simulation thread
	void OnEndFrame() override;

private:
	stdc::time_point<stdclock> m_SimFrameEnd;
};

