#include "Test.hpp"
#include "Errore.hpp"
#include "Utility.hpp"

#include "OpenMP/SPHSimulation.hpp"
#include "Serial/SPHSimulation.hpp"

#include "Visualizers/FileExporter.hpp"
#ifndef DISABLE_UI
#include "Visualizers/ImGuiViewer.hpp"
#endif

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"
#include "Initializers/PipeInitializer.hpp"
#include "Initializers/PascalInitializer.hpp"

#include "OpenMP/Neighbors/SpatialHashing.hpp"
#include "OpenMP/Neighbors/MortonSorting.hpp"

#include "Serial/Neighbors/SpatialHashing.hpp"
#include "Serial/Neighbors/MortonSorting.hpp"

#ifndef DISABLE_UI
void Test_Correctness_Check()
{
	CompareSeries<2>("Correctness/AoS_OMP_Morton", "Correctness/AoS_Serial_Morton");
	CompareSeries<2>("Correctness/AoS_OMP_Morton", "Correctness/SoA_OMP_Morton");
	CompareSeries<2>("Correctness/AoS_OMP_Morton", "Correctness/Hybrid_OMP_Morton");
}
void Test_Correctness_Generate()
{
	if (!fs::exists("Correctness/"))
		fs::create_directories("Correctness");

	base::SPHParams params;
	params.FinalTime = params.TimeStep * 100;

	GENERATE_TEST_CASE_COMPARE(2, openmp, MortonSorting, ParticleAoS, "Correctness/AoS_OMP_Morton", BoxInitializer, "AoS_OMP_Morton", params);
	GENERATE_TEST_CASE_COMPARE(2, serial, MortonSorting, ParticleAoS, "Correctness/AoS_Serial_Morton", BoxInitializer, "AoS_Serial_Morton", params);
	GENERATE_TEST_CASE_COMPARE(2, openmp, SpatialHashing, ParticleAoS, "Correctness/AoS_OMP_Hash", BoxInitializer, "AoS_OMP_Hash", params);
	GENERATE_TEST_CASE_COMPARE(2, serial, SpatialHashing, ParticleAoS, "Correctness/AoS_Serial_Hash", BoxInitializer, "AoS_Serial_Hash", params);

	GENERATE_TEST_CASE_COMPARE(2, openmp, MortonSorting, ParticleHybrid, "Correctness/Hybrid_OMP_Morton", BoxInitializer, "Hybrid_OMP_Morton", params);
	GENERATE_TEST_CASE_COMPARE(2, serial, MortonSorting, ParticleHybrid, "Correctness/Hybrid_Serial_Morton", BoxInitializer, "Hybrid_Serial_Morton", params);
	GENERATE_TEST_CASE_COMPARE(2, openmp, SpatialHashing, ParticleHybrid, "Correctness/Hybrid_OMP_Hash", BoxInitializer, "Hybrid_OMP_Hash", params);
	GENERATE_TEST_CASE_COMPARE(2, serial, SpatialHashing, ParticleHybrid, "Correctness/Hybrid_Serial_Hash", BoxInitializer, "Hybrid_Serial_Hash", params);
	
	GENERATE_TEST_CASE_COMPARE(2, openmp, MortonSorting, ParticleSoA, "Correctness/SoA_OMP_Morton", BoxInitializer, "SoA_OMP_Morton", params);
	GENERATE_TEST_CASE_COMPARE(2, serial, MortonSorting, ParticleSoA, "Correctness/SoA_Serial_Morton", BoxInitializer, "SoA_Serial_Morton", params);
	GENERATE_TEST_CASE_COMPARE(2, openmp, SpatialHashing, ParticleSoA, "Correctness/SoA_OMP_Hash", BoxInitializer, "SoA_OMP_Hash", params);
	GENERATE_TEST_CASE_COMPARE(2, serial, SpatialHashing, ParticleSoA, "Correctness/SoA_Serial_Hash", BoxInitializer, "SoA_Serial_Hash", params);
}

void Test_2D_OpenMP_Morton()
{
	GENERATE_TEST_CASE_VIEW(2, openmp, MortonSorting, ParticleAoS, BoxInitializer, "2D_OMP_AoS_Morton", base::SPHParams{});
}
#endif
void Test_2D_OpenMP_Hashing()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	openmp::SpatialHashing<Size, Particles> nf;

	openmp::SPHSimulation<Size, Particles> sph(&nf);
	sph.SetName("AoS-OMP-Hashing");

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

#ifndef DISABLE_UI
	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
#endif
	sph.Start();
}

void Test_2D_Serial_Morton()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	serial::MortonSorting<Size, Particles> nf({ 0.0f, 0.0f }, { 1.0f, 1.0f });

	serial::SPHSimulation<Size, Particles> sph(&nf);
	sph.SetName("AoS-Serial-Morton");

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

#ifndef DISABLE_UI
	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
#endif
	sph.Start();
}
void Test_2D_Serial_Hashing()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	serial::SpatialHashing<Size, Particles> nf;

	serial::SPHSimulation<Size, Particles> sph(&nf);
	sph.SetName("AoS-Serial-Hashing");

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

#ifndef DISABLE_UI
	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
#endif
	sph.Start();
}
void Test_2D_Pipe_OpenMP_Morton()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	openmp::MortonSorting<Size, Particles> nf({ 0.0f, 0.0f }, { 1.0f, 1.0f });

	openmp::SPHSimulation<Size, Particles> sph(&nf);
	sph.SetName("AoS-Pipe-OpenMP-Morton");

	PipeInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

#ifndef DISABLE_UI
	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
#endif
	sph.Start();

}
void Test_2D_Pascal_OpenMP_Morton()
{
	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	openmp::MortonSorting<Size, Particles> nf({ 0.0f, 0.0f }, { 1.0f, 1.0f });

	openmp::SPHSimulation<Size, Particles> sph(&nf);
	sph.SetName("AoS-Pascal-OpenMP-Morton");

	PascalInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);
#ifndef DISABLE_UI
	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
#endif
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
