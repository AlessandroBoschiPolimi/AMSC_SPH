
#include "Utility.hpp"
#include "SPHSimulation.hpp"
#include "XYZExporter.hpp"
#include "ImGuiViewer.hpp"


#define SIZE 3

// Main code
int main(int, char**)
{
	SPHSimulation<SIZE> sph;


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