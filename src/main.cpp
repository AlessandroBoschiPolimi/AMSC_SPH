
#include "Utility.hpp"
#include "OpenMP/SPHSimulation.hpp"
#include "Serial/SPHSimulation.hpp"

#include "Visualizers/FileExporter.hpp"
#include "Visualizers/ImGuiViewer.hpp"

#include "OpenMP/Neighbors/SpatialHashing.hpp"
#include "OpenMP/Neighbors/MortonSorting.hpp"

#include "Serial/Neighbors/SpatialHashing.hpp"
#include "Serial/Neighbors/MortonSorting.hpp"

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"

#include <cstring>

#include "CUDA/Test.cuh"


#define SIZE 2

void ParseArgs(int argc, char** argv);

void Test_2D_OpenMP_Morton();
void Test_2D_OpenMP_Hashing();
void Test_2D_Serial_Morton();
void Test_2D_Serial_Hashing();


// Main code
int main(int argc, char** argv)
{
#ifdef HAS_CUDA
	Test();
#else
	std::cout << "NO CUDA\n";
#endif

	ParseArgs(argc, argv);
	
	//omp_set_num_threads(8);
	std::cout << "OpenMP Max Threads: " << omp_get_max_threads() << '\n';

	Test_2D_OpenMP_Morton();

	return 0;
}


void Test_2D_OpenMP_Morton()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	// TODO: define first the initializer, and query the domain min and max from it
	openmp::MortonSorting<Size, Particles> nf({ 0.0f, 0.0f }, { 1.0f, 1.0f });

	openmp::SPHSimulation<Size, Particles> sph(&nf);

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();

	sph.Start();
}
void Test_2D_OpenMP_Hashing()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	openmp::SpatialHashing<Size, Particles> nf;

	openmp::SPHSimulation<Size, Particles> sph(&nf);

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();

	sph.Start();
}

void Test_2D_Serial_Morton()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	// TODO: define first the initializer, and query the domain min and max from it
	serial::MortonSorting<Size, Particles> nf({ 0.0f, 0.0f }, { 1.0f, 1.0f });

	serial::SPHSimulation<Size, Particles> sph(&nf);

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();

	sph.Start();
}
void Test_2D_Serial_Hashing()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	serial::SpatialHashing<Size, Particles> nf;

	serial::SPHSimulation<Size, Particles> sph(&nf);

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();

	sph.Start();
}

//void Test_3D_OpenMP_Morton()
//{
//	constexpr size_t Size = 3;
//	using Particles = ParticleAoS<Size>;
//
//	openmp::MortonSorting<Size, Particles> nf({ 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f });
//
//	openmp::SPHSimulation<Size, Particles> sph(&nf);
//
//	TrayInitializer<Size, Particles> initializer;
//	sph.InitializeFluid(&initializer);
//
//	FileExporter<Particles> exporter;
//	sph.AddObserver(&exporter);
//
//	sph.Start();
//}

void ParseArgs(int argc, char** argv)
{
	for (int i = 1; i < argc; true)
	{
		char* line = argv[i++];
		if (std::strcmp(line, "--omp-threads") == 0 && i < argc)
		{
			int count = std::stoi(argv[i++]);
			omp_set_num_threads(count);
		}
	}
}