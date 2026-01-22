
#include "Utility.hpp"
#include "SPHSimulation.hpp"
#include "XYZExporter.hpp"
#include "ImGuiViewer.hpp"



// Main code
int main(int, char**)
{
	SPHSimulation sph;

	ImGuiViewer imguiViewer;
	//XYZExporter exporter;

	sph.AddObserver(&imguiViewer);
	//sph.AddObserver(&exporter);

	imguiViewer.Start();
	sph.Start();

	return 0;
}