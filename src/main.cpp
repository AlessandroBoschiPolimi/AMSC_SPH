#include "Utility.hpp"
#include "Test.hpp"

#include "OpenMP/SPHSimulation.hpp"
#include "CUDA/SPHSimulation.hpp"
#include "Serial/SPHSimulation.hpp"

#include "Visualizers/FileExporter.hpp"
#include "Visualizers/ImGuiViewer.hpp"
#include "Visualizers/CudaImGuiViewer.hpp"

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"

#include "OpenMP/Neighbors/SpatialHashing.hpp"
#include "OpenMP/Neighbors/MortonSorting.hpp"

#include "Serial/Neighbors/SpatialHashing.hpp"
#include "Serial/Neighbors/MortonSorting.hpp"

#include <cstring>


void ParseArgs(int argc, char** argv);

void TestCuda();

int main(int argc, char** argv)
{
	ParseArgs(argc, argv);
	
	// omp_set_num_threads(8);
	std::cout << "OpenMP Max Threads: " << omp_get_max_threads() << '\n';

	// Test_Correctness_Generate();
	// Test_Correctness_Check();

	// base::SPHParams params;
	// params.FinalTime = params.TimeStep * 100;

	// GENERATE_TEST_CASE_PARAVIEW(2, openmp, MortonSorting, ParticleAoS, "Correctness/AoS_OMP_Morton", BoxInitializer, "AoS_OMP_Morton", params);
	// GENERATE_TEST_CASE_PARAVIEW(2, openmp, MortonSorting, ParticleSoA, "Correctness/SoA_OMP_Morton", BoxInitializer, "SoA_OMP_Morton", params);

	// GENERATE_TEST_CASE_VIEW(2, openmp, MortonSorting, ParticleAoS, BoxInitializer, "2D_OMP_AoS_Morton", base::SPHParams{});

	TestCuda();

	return 0;
}

void TestCuda()
{
#ifdef HAS_CUDA
	using Particles = ParticleSoA<2>;
	BoxInitializer<2, Particles> init;

	cudasph::SPHSimulation<2> sph;
	sph.SetName("CUDA");
	sph.SetParams(base::SPHParams{});
	sph.InitializeFluid(&init);

	CudaImGuiViewer viewer;
	sph.AddObserver(&viewer);
	viewer.Start();

	sph.Start();
#else
	std::cout << "NO CUDA\n";
#endif
}


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