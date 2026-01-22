
#include "Utility.hpp"
#include "SPHSimulation.hpp"
#include "XYZExporter.hpp"
#include "ImGuiViewer.hpp"



// Main code
int main(int, char**)
{
	SPHSimulation sph;

	ImGuiViewer imguiViewer;
	// XYZExporter exporter;

	sph.AddObserver(&imguiViewer);
	// sph.addObserver(&exporter);

	imguiViewer.Start();
	sph.Start();

	std::cout << "End\n";

	return 0;
}