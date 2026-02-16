
#include "Utility.hpp"
#include "OpenMP/SPHSimulation.hpp"
#include "Serial/SPHSimulation.hpp"

#include "Visualizers/XYZExporter.hpp"
#include "Visualizers/ImGuiViewer.hpp"

#include "OpenMP/Neighbors/SpatialHashing.hpp"
#include "OpenMP/Neighbors/MortonSorting.hpp"

#include "Serial/Neighbors/SpatialHashing.hpp"
#include "Serial/Neighbors/MortonSorting.hpp"

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"

#include <cstring>



#define SIZE 2

void ParseArgs(int argc, char** argv);


// Main code
int main(int argc, char** argv)
{
	ParseArgs(argc, argv);
	
	std::cout << "OpenMP Max Threads: " << omp_get_max_threads() << '\n';

#if SIZE == 2
	// TODO: define first the initializer, and query the domain min and max from it
	openmp::MortonSorting<SIZE> nf({ 0.0f, 0.0f }, { 1.0f, 1.0f });
#else
	openmp::MortonSorting<SIZE> nf({ 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f });
#endif

	openmp::SPHSimulation<SIZE> sph(&nf);


	//std::cout << fs::current_path() << '\n';

#if SIZE == 2
	BoxInitializer<SIZE> initializer;
	sph.InitializeFluid(&initializer);

	ImGuiViewer imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
#else
	TrayInitializer<SIZE> initializer;
	sph.InitializeFluid(&initializer);

	XYZExporter exporter;
	sph.AddObserver(&exporter);
#endif

	sph.Start();

	return 0;
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