#pragma once
#include "Observer.hpp"
#include <fstream>
#include <format>

void writeXYZ(const int frame, const std::vector<Particle<3>>& particles);
void writeVTU(const int frame, const std::vector<Particle<3>>& particles);
void writeVTUBinary(const int frame, const std::vector<Particle<3>>& particles);

class XYZExporter : public Observer<3> {
public:
	~XYZExporter() override = default;

	/// Executes on the simulation thread
	void OnEndFrame() override {
		if (m_Sim == nullptr)
			return;

		std::cout << "Frame: " << m_Sim->GetFrame() << '\n';

		auto now = stdclock::now();
		std::cout << "Simulate ";
		print_time(now - m_SimFrameEnd);
		std::cout << ' ';
		
		auto writeStart = stdclock::now();
		writeVTUBinary(m_Sim->GetFrame(), m_Sim->GetParticles());
		auto writeEnd = stdclock::now();
		std::cout << "Write ";
		print_time(writeEnd - writeStart);
		std::cout << '\n';

		m_SimFrameEnd = stdclock::now();
	}

private:
	stdc::time_point<stdclock> m_SimFrameEnd;
};

void writeXYZ(const int frame, const std::vector<Particle<3>>& particles)
{
	std::ofstream out(std::format("output-{}.xyz", frame));
	for (auto& p : particles)
		out << p.Position.x << " " << p.Position.y << " " << p.Position.z << "\n";
}

void writeVTU(const int frame, const std::vector<Particle<3>>& particles) {
	std::ofstream out(std::format("output-{}.vtu", frame));
	const size_t N = particles.size();

	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "  <UnstructuredGrid>\n";
	out << "    <Piece NumberOfPoints=\"" << N
		<< "\" NumberOfCells=\"" << N << "\">\n";

	// --- Points ---
	out << "      <Points>\n";
	out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (const auto& p : particles)
		out << p.Position.x << " " << p.Position.y << " " << p.Position.z << "\n";
	out << "        </DataArray>\n";
	out << "      </Points>\n";

	// --- Point Data (scalars, vectors) ---
	out << "      <PointData Scalars=\"density\">\n";

	// Density
	out << "        <DataArray type=\"Float32\" Name=\"density\" format=\"ascii\">\n";
	for (const auto& p : particles)
		out << p.Density << "\n";
	out << "        </DataArray>\n";

	// Pressure
	out << "        <DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">\n";
	for (const auto& p : particles)
		out << p.Pressure << "\n";
	out << "        </DataArray>\n";

	// Velocity (vector)
	out << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (const auto& p : particles)
		out << p.Velocity.x << " " << p.Velocity.y << " " << p.Velocity.z << "\n";
	out << "        </DataArray>\n";

	out << "      </PointData>\n";

	// --- Cells ---
	out << "      <Cells>\n";

	// Connectivity
	out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	for (size_t i = 0; i < N; ++i)
		out << i << "\n";
	out << "        </DataArray>\n";

	// Offsets
	out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for (size_t i = 1; i <= N; ++i)
		out << i << "\n";
	out << "        </DataArray>\n";

	// Types (1 = VTK_VERTEX)
	out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for (size_t i = 0; i < N; ++i)
		out << "1\n";
	out << "        </DataArray>\n";

	out << "      </Cells>\n";

	out << "    </Piece>\n";
	out << "  </UnstructuredGrid>\n";
	out << "</VTKFile>\n";
}

template <typename T>
void writeRaw(std::ofstream& out, const T& value) {
	out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

void writeVTUBinary(const int frame, const std::vector<Particle<3>>& particles) {
	std::ofstream out(std::format("output-{}.vtu", frame), std::ios::binary);
	if (!out) return;

	const uint32_t N = static_cast<uint32_t>(particles.size());

	// --- Compute byte sizes ---
	uint32_t pointsSize = N * 3 * sizeof(float);
	uint32_t densitySize = N * sizeof(float);
	uint32_t pressureSize = N * sizeof(float);
	uint32_t velocitySize = N * 3 * sizeof(float);
	uint32_t connSize = N * sizeof(int32_t);
	uint32_t offsetsSize = N * sizeof(int32_t);
	uint32_t typesSize = N * sizeof(uint8_t);

	// --- Compute offsets (relative to start of appended data, after '_') ---
	uint32_t offset = 0;

	uint32_t pointsOffset = offset; offset += sizeof(uint32_t) + pointsSize;
	uint32_t densityOffset = offset; offset += sizeof(uint32_t) + densitySize;
	uint32_t pressureOffset = offset; offset += sizeof(uint32_t) + pressureSize;
	uint32_t velocityOffset = offset; offset += sizeof(uint32_t) + velocitySize;
	uint32_t connOffset = offset; offset += sizeof(uint32_t) + connSize;
	uint32_t offsetsOffset = offset; offset += sizeof(uint32_t) + offsetsSize;
	uint32_t typesOffset = offset; offset += sizeof(uint32_t) + typesSize;

	// --- XML Header ---
	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
		<< "byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	out << "  <UnstructuredGrid>\n";
	out << "    <Piece NumberOfPoints=\"" << N
		<< "\" NumberOfCells=\"" << N << "\">\n";

	// Points
	out << "      <Points>\n";
	out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" "
		<< "format=\"appended\" offset=\"" << pointsOffset << "\"/>\n";
	out << "      </Points>\n";

	// Point data
	out << "      <PointData Scalars=\"density\">\n";

	out << "        <DataArray type=\"Float32\" Name=\"density\" "
		<< "format=\"appended\" offset=\"" << densityOffset << "\"/>\n";

	out << "        <DataArray type=\"Float32\" Name=\"pressure\" "
		<< "format=\"appended\" offset=\"" << pressureOffset << "\"/>\n";

	out << "        <DataArray type=\"Float32\" Name=\"velocity\" "
		<< "NumberOfComponents=\"3\" format=\"appended\" "
		<< "offset=\"" << velocityOffset << "\"/>\n";

	out << "      </PointData>\n";

	// Cells
	out << "      <Cells>\n";

	out << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
		<< "format=\"appended\" offset=\"" << connOffset << "\"/>\n";

	out << "        <DataArray type=\"Int32\" Name=\"offsets\" "
		<< "format=\"appended\" offset=\"" << offsetsOffset << "\"/>\n";

	out << "        <DataArray type=\"UInt8\" Name=\"types\" "
		<< "format=\"appended\" offset=\"" << typesOffset << "\"/>\n";

	out << "      </Cells>\n";

	out << "    </Piece>\n";
	out << "  </UnstructuredGrid>\n";

	// --- Appended Data ---
	out << "  <AppendedData encoding=\"raw\">\n_";

	// Points
	writeRaw(out, pointsSize);
	for (const auto& p : particles) {
		writeRaw(out, p.Position.x);
		writeRaw(out, p.Position.y);
		writeRaw(out, p.Position.z);
	}

	// Density
	writeRaw(out, densitySize);
	for (const auto& p : particles)
		writeRaw(out, p.Density);

	// Pressure
	writeRaw(out, pressureSize);
	for (const auto& p : particles)
		writeRaw(out, p.Pressure);

	// Velocity
	writeRaw(out, velocitySize);
	for (const auto& p : particles) {
		writeRaw(out, p.Velocity.x);
		writeRaw(out, p.Velocity.y);
		writeRaw(out, p.Velocity.z);
	}

	// Connectivity
	writeRaw(out, connSize);
	for (uint32_t i = 0; i < N; ++i)
		writeRaw(out, static_cast<int32_t>(i));

	// Offsets
	writeRaw(out, offsetsSize);
	for (uint32_t i = 1; i <= N; ++i)
		writeRaw(out, static_cast<int32_t>(i));

	// Types (VTK_VERTEX = 1)
	writeRaw(out, typesSize);
	for (uint32_t i = 0; i < N; ++i)
		writeRaw(out, static_cast<uint8_t>(1));

	out << "\n  </AppendedData>\n";
	out << "</VTKFile>\n";
}