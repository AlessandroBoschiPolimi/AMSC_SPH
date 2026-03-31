#include "Utility.hpp"

#include "Particle.hpp"


struct Stats
{
	double max = 0.0;
	double sum = 0.0;
	size_t count = 0;

	void Add(double v)
	{
		max = std::max(max, v);
		sum += v;
		++count;
	}

	double Avg() const
	{
		return count ? sum / count : 0.0;
	}
};

struct FieldStats
{
	Stats position;
	Stats velocity;
	Stats a_grav;
	Stats a_press;
	Stats a_visc;
	Stats density;
	Stats pressure;
};

template <size_t D>
Particle<D> ReadParticle(std::ifstream& in);

template <>
inline Particle<2> ReadParticle(std::ifstream& in)
{
	Particle<2> p;
	in >> p.Position.x >> p.Position.y;
	in >> p.Velocity.x >> p.Velocity.y;
	in >> p.A_grav.x >> p.A_grav.y;
	in >> p.A_press.x >> p.A_press.y;
	in >> p.A_visc.x >> p.A_visc.y;
	in >> p.Density;
	in >> p.Pressure;
	return p;
}
template <>
inline Particle<3> ReadParticle(std::ifstream& in)
{
	Particle<3> p;
	in >> p.Position.x >> p.Position.y >> p.Position.z;
	in >> p.Velocity.x >> p.Velocity.y >> p.Velocity.z;
	in >> p.A_grav.x >> p.A_grav.y >> p.A_grav.z;
	in >> p.A_press.x >> p.A_press.y >> p.A_press.z;
	in >> p.A_visc.x >> p.A_visc.y >> p.A_visc.z;
	in >> p.Density;
	in >> p.Pressure;
	return p;
}

template <size_t D>
std::vector<Particle<D>> ReadFile(const std::string& filename)
{
	std::ifstream in(filename);

	size_t size;
	in >> size;

	std::vector<Particle<D>> particles;
	particles.resize(size);

	for (size_t i = 0; i < size; i++)
		particles[i] = ReadParticle<D>(in);

	in.close();

	return particles;
}

template <size_t D>
void CompareFrame(const std::string& fileA, const std::string& fileB, FieldStats& globalStats)
{
	std::vector<Particle<D>> vecA = ReadFile<D>(fileA);
	std::vector<Particle<D>> vecB = ReadFile<D>(fileB);

	if (vecA.size() != vecB.size())
		throw std::runtime_error("Uncompatible file sizes");

	for (size_t i = 0; i < vecA.size(); ++i)
	{
		Particle<D> pA = vecA[i], pB = vecB[i];

		globalStats.position.Add(Norm(pA.Position - pB.Position));
		globalStats.velocity.Add(Norm(pA.Velocity - pB.Velocity));
		globalStats.a_grav.Add(Norm(pA.A_grav - pB.A_grav));
		globalStats.a_visc.Add(Norm(pA.A_visc - pB.A_visc));
		globalStats.a_press.Add(Norm(pA.A_press - pB.A_press));
		globalStats.density.Add(pA.Density - pB.Density);
		globalStats.pressure.Add(pA.Pressure - pB.Pressure);
	}
}

template <size_t D>
void CompareSeries(const std::string& baseA, const std::string& baseB)
{
	FieldStats stats;

	u64 frame = 1;
	while (true)
	{
		std::string fileA, fileB;

		{
			std::ostringstream oss;
			oss << baseA << '-' << frame << ".mycoolextension";
			fileA = oss.str();
			oss.str("");
			oss << baseB << '-' << frame << ".mycoolextension";
			fileB = oss.str();
		}

		if (!std::filesystem::exists(fileA) || !std::filesystem::exists(fileB))
			break;

		std::cout << "Comparing frame " << frame << '\n';

		CompareFrame<D>(fileA, fileB, stats);
		++frame;
	}

	std::cout << "\n==== Global Results Across " << frame << " Frames ====\n";

	std::cout << "Position:    max = " << stats.position.max << ", avg = " << stats.position.Avg() << "\n";

	std::cout << "Velocity:    max = " << stats.velocity.max << ", avg = " << stats.velocity.Avg() << "\n";

	std::cout << "Acceleration: max = " << stats.a_grav.max << ", avg = " << stats.a_grav.Avg() << "\n";
	std::cout << "Acceleration: max = " << stats.a_press.max << ", avg = " << stats.a_press.Avg() << "\n";
	std::cout << "Acceleration: max = " << stats.a_visc.max << ", avg = " << stats.a_visc.Avg() << "\n";

	std::cout << "Density:     max = " << stats.density.max << ", avg = " << stats.density.Avg() << "\n";

	std::cout << "Pressure:    max = " << stats.pressure.max << ", avg = " << stats.pressure.Avg() << "\n";
	
	std::cout << '\n';
}