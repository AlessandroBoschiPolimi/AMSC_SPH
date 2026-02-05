
#include "Utility.hpp"
#include "SPHSimulation.hpp"

#include "Visualizers/XYZExporter.hpp"
#include "Visualizers/ImGuiViewer.hpp"

#include "Neighbors/SpatialHashing.hpp"

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"


#define SIZE 2

// Main code
int main(int, char**)
{
	SpatialHashing<SIZE> nf;

	SPHSimulation<SIZE> sph(&nf);


	std::cout << fs::current_path() << '\n';

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