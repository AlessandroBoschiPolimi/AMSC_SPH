
#include "Utility.hpp"
#include "SPHSimulation.hpp"
#include "XYZExporter.hpp"
#include "ImGuiViewer.hpp"


#define SIZE 2

// Main code
int main(int, char**)
{
	SPHSimulation<SIZE> sph;
	BoxInitializer<SIZE> initializer;

	sph.InitializeFluid(&initializer);

	std::cout << fs::current_path();

#if SIZE == 2
	ImGuiViewer imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
#else
	XYZExporter exporter;
	sph.AddObserver(&exporter);
#endif

	sph.Start();

	return 0;
}