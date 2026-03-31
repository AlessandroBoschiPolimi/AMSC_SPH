#include <mpi.h>

#include "Utility.hpp"
#include "Test.hpp"

#include "OpenMP/SPHSimulation.hpp"
#include "Serial/SPHSimulation.hpp"
#include "MPI/SPHSimulation.hpp"

#include "Visualizers/FileExporter.hpp"
#ifndef DISABLE_UI
#include "Visualizers/ImGuiViewer.hpp"
#endif

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"

#include "OpenMP/Neighbors/SpatialHashing.hpp"
#include "OpenMP/Neighbors/MortonSorting.hpp"

#include "Serial/Neighbors/SpatialHashing.hpp"
#include "Serial/Neighbors/MortonSorting.hpp"

#include "MPI/Neighbors/SpatialHashing.hpp"
#include "MPI/Neighbors/MortonSorting.hpp"


#include <cstring>

#define SIZE 2

void ParseArgs(int argc, char** argv);
void ParseExamplesArgs(int argc, char** argv);


// Main code
int main(int argc, char** argv)
{
	ParseArgs(argc, argv);
	
	std::cout << "OpenMP Max Threads: " << omp_get_max_threads() << '\n';
	
	ParseExamplesArgs(argc, argv);
	
	// omp_set_num_threads(8);
	Test_2D_OpenMP_Hashing();
	// Test_Correctness_Generate();
	// Test_Correctness_Check();

	// base::SPHParams params;
	// params.FinalTime = params.TimeStep * 100;

	// GENERATE_TEST_CASE_PARAVIEW(2, openmp, MortonSorting, ParticleAoS, "Correctness/AoS_OMP_Morton", BoxInitializer, "AoS_OMP_Morton", params);
	// GENERATE_TEST_CASE_PARAVIEW(2, openmp, MortonSorting, ParticleSoA, "Correctness/SoA_OMP_Morton", BoxInitializer, "SoA_OMP_Morton", params);

	// GENERATE_TEST_CASE_VIEW(2, openmp, MortonSorting, ParticleAoS, BoxInitializer, "2D_OMP_AoS_Morton", base::SPHParams{});

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

void ParseExamplesArgs(int argc, char** argv)
{
	int size = 2, par = 0, init = 0, nf = 0, layout = 0;
	bool file = false, ui = false;
	std::string filename;
	ExportFormat fileformat = ExportFormat::VTU_01;

	bool run_example = false;

	for (int i = 1; i < argc; true)
	{
		char* line = argv[i++];
		if (std::strcmp(line, "-example") == 0)
			run_example = true;
		else if (std::strcmp(line, "-d") == 0)
		{
			if (i >= argc) {
				std::cerr << "Provide either '2' or '3' after flag '-d'";
				exit(1);
			}

			int val = std::stoi(argv[i++]);
			size = val;
			if (size != 2 && size != 3) {
				std::cerr << "Provide either '2' or '3' after flag '-d'";
				exit(1);
			}
		}
		else if (std::strcmp(line, "-par") == 0)
		{
			if (i >= argc)
			{
				std::cerr << "Provide either 'serial', 'openmp' or 'mpi' after flag '-par'";
				exit(1);
			}

			line = argv[i++];
			if (std::strcmp(line, "serial") == 0)
				par = 0;
			else if (std::strcmp(line, "openmp") == 0)
				par = 1;
			else if (std::strcmp(line, "mpi") == 0)
				par = 2;
			else {
				std::cerr << "Provide either 'serial', 'openmp' or 'mpi' after flag '-par'";
				exit(1);
			}
		}
		else if (std::strcmp(line, "-layout") == 0)
		{
			if (i >= argc)
			{
				std::cerr << "Provide either 'AoS', 'SoA' or 'hybrid' after flag '-layout'";
				exit(1);
			}

			line = argv[i++];
			if (std::strcmp(line, "AoS") == 0)
				layout = 0;
			else if (std::strcmp(line, "SoA") == 0)
				layout = 1;
			else if (std::strcmp(line, "hybrid") == 0)
				layout = 2;
			else {
				std::cerr << "Provide either 'AoS', 'SoA' or 'hybrid' after flag '-layout'";
				exit(1);
			}
		}
		else if (std::strcmp(line, "-init") == 0)
		{
			if (i >= argc)
			{
				std::cerr << "Provide either 'box', 'tray' or 'TODO' after flag '-init'";
				exit(1);
			}

			line = argv[i++];
			if (std::strcmp(line, "box") == 0)
				init = 0;
			else if (std::strcmp(line, "tray") == 0)
				init = 1;
			else {
				std::cerr << "Provide either 'box', 'tray' or 'TODO' after flag '-init'";
				exit(1);
			}
		}
		else if (std::strcmp(line, "-nf") == 0)
		{
			if (i >= argc)
			{
				std::cerr << "Provide either 'hash', or 'morton' after flag '-nf'";
				exit(1);
			}

			line = argv[i++];
			if (std::strcmp(line, "hash") == 0)
				nf = 0;
			else if (std::strcmp(line, "morton") == 0)
				nf = 1;
			else {
				std::cerr << "Provide either 'hash', or 'morton' after flag '-nf'";
				exit(1);
			}
		}
		else if (std::strcmp(line, "-file") == 0)
		{
			if (i + 1 >= argc)
			{
				std::cerr << "Provide filename and either 'xyz', 'vtu' or 'bvtu' after flag '-file'";
				exit(1);
			}

			filename = std::string(argv[i++]);

			line = argv[i++];
			if (std::strcmp(line, "xyz") == 0)
				fileformat = ExportFormat::XYZ;
			else if (std::strcmp(line, "vtu") == 0)
				fileformat = ExportFormat::VTU;
			else if (std::strcmp(line, "bvtu") == 0)
				fileformat = ExportFormat::VTU_01;
			else {
				std::cerr << "Provide filename and either 'xyz', 'vtu' or 'bvtu' after flag '-file'";
				exit(1);
			}

			file = true;
		}
		else if (std::strcmp(line, "-ui") == 0)
		{
			ui = true;
		}
	}


	if (!run_example)
		return;

	if (size == 2)
	{
		static constexpr int Size = 2;

		if (par == 2)
		{
			MPI_Init(&argc, &argv);
			MPI_Comm mpi_comm = MPI_COMM_WORLD;
			//Get the number of processes .
			int mpi_size;
			MPI_Comm_size(mpi_comm, &mpi_size);
			//Get the rank of current process .
			int mpi_rank;
			MPI_Comm_rank(mpi_comm, &mpi_rank);

			#ifndef DISABLE_UI
			ImGuiViewer<ParticleAoS<Size>> imguiViewer;
			#endif
			FileExporter<Size, ParticleAoS<Size>> exporter;

			SimInitializer<Size, ParticleAoS<Size>>* initp;

			switch (init)
			{
			default: initp = new BoxInitializer<Size, ParticleAoS<Size>>(); break;
			}

			mpi::NeighborFinder<Size>* nfp;

			switch (nf)
			{
			case 0: nfp = new mpi::SpatialHashing<Size>(); break;
			default: nfp = new mpi::MortonSorting<Size>(initp->GetDomain().first, initp->GetDomain().second); break;
			}

			mpi::SPHSimulation<Size> sim(nfp);

			if (nf == 0)
				sim.SetName("MPI AoS Hashing");
			else
				sim.SetName("MPI AoS Morton");

			sim.SetRank(mpi_rank, mpi_comm, mpi_size);
			if (mpi_rank == 0)
			{
				sim.InitializeFluid(initp);

				#ifndef DISABLE_UI
				if (ui)
				{
					sim.AddObserver(&imguiViewer);
					imguiViewer.Start();
				}
				#endif
				if (file)
				{
					exporter.SetBaseName(filename);
					exporter.SetFormat(fileformat);
					sim.AddObserver(&exporter);
				}
			}

			sim.Start();

			delete nfp;
			delete initp;

			// Finalize MPI.
			MPI_Finalize();
		}

		else
		{
			if (layout == 0)
			{
				using Particles = ParticleAoS<Size>;
				#ifndef DISABLE_UI
				ImGuiViewer<Particles> imguiViewer;
				#endif
				FileExporter<Size, Particles> exporter;

				SimInitializer<Size, Particles>* initp;

				switch (init)
				{
				default: initp = new BoxInitializer<Size, Particles>(); break;
				}

				if (par == 0)
				{
					serial::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new serial::SpatialHashing<Size, Particles>(); break;
					default: nfp = new serial::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					serial::SPHSimulation<Size, Particles> sim(nfp);
					
					if (nf == 0)
						sim.SetName("Serial AoS Hashing");
					else
						sim.SetName("Serial AoS Morton");
					
					sim.InitializeFluid(initp);

					if (ui)
					{
						#ifndef DISABLE_UI
						sim.AddObserver(&imguiViewer);
						imguiViewer.Start();
						#endif
					}
					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				else if (par == 1)
				{
					openmp::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new openmp::SpatialHashing<Size, Particles>(); break;
					default: nfp = new openmp::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					openmp::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("OMP AoS Hashing");
					else
						sim.SetName("OMP AoS Morton");

					sim.InitializeFluid(initp);

					if (ui)
					{
						#ifndef DISABLE_UI
						sim.AddObserver(&imguiViewer);
						imguiViewer.Start();
						#endif
					}
					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				delete initp;
			}
			
			else if (layout == 1)
			{
				using Particles = ParticleSoA<Size>;
				#ifndef DISABLE_UI
				ImGuiViewer<Particles> imguiViewer;
				#endif
				FileExporter<Size, Particles> exporter;

				SimInitializer<Size, Particles>* initp;

				switch (init)
				{
				default: initp = new BoxInitializer<Size, Particles>(); break;
				}

				if (par == 0)
				{
					serial::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new serial::SpatialHashing<Size, Particles>(); break;
					default: nfp = new serial::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					serial::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("Serial SoA Hashing");
					else
						sim.SetName("Serial SoA Morton");

					sim.InitializeFluid(initp);

					if (ui)
					{
						#ifndef DISABLE_UI
						sim.AddObserver(&imguiViewer);
						imguiViewer.Start();
						#endif
					}
					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				else if (par == 1)
				{
					openmp::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new openmp::SpatialHashing<Size, Particles>(); break;
					default: nfp = new openmp::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					openmp::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("OMP SoA Hashing");
					else
						sim.SetName("OMP SoA Morton");

					sim.InitializeFluid(initp);

					if (ui)
					{
						#ifndef DISABLE_UI
						sim.AddObserver(&imguiViewer);
						imguiViewer.Start();
						#endif
					}
					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				delete initp;
			}

			else if (layout == 2)
			{
				using Particles = ParticleHybrid<Size>;
				#ifndef DISABLE_UI
				ImGuiViewer<Particles> imguiViewer;
				#endif
				FileExporter<Size, Particles> exporter;

				SimInitializer<Size, Particles>* initp;

				switch (init)
				{
				default: initp = new BoxInitializer<Size, Particles>(); break;
				}

				if (par == 0)
				{
					serial::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new serial::SpatialHashing<Size, Particles>(); break;
					default: nfp = new serial::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					serial::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("Serial Hybrid Hashing");
					else
						sim.SetName("Serial Hybrid Morton");

					sim.InitializeFluid(initp);

					if (ui)
					{
						#ifndef DISABLE_UI
						sim.AddObserver(&imguiViewer);
						imguiViewer.Start();
						#endif
					}
					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				else if (par == 1)
				{
					openmp::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new openmp::SpatialHashing<Size, Particles>(); break;
					default: nfp = new openmp::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					openmp::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("OMP Hybrid Hashing");
					else
						sim.SetName("OMP Hybrid Morton");

					sim.InitializeFluid(initp);

					if (ui)
					{
						#ifndef DISABLE_UI
						sim.AddObserver(&imguiViewer);
						imguiViewer.Start();
						#endif
					}
					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				delete initp;
			}
		}
	}
	else
	{
		static constexpr int Size = 3;

		if (par == 2)
		{
			MPI_Init(&argc, &argv);
			MPI_Comm mpi_comm = MPI_COMM_WORLD;
			//Get the number of processes .
			int mpi_size;
			MPI_Comm_size(mpi_comm, &mpi_size);
			//Get the rank of current process .
			int mpi_rank;
			MPI_Comm_rank(mpi_comm, &mpi_rank);

			FileExporter<Size, ParticleAoS<Size>> exporter;

			SimInitializer<Size, ParticleAoS<Size>>* initp;

			switch (init)
			{
			default: initp = new TrayInitializer<Size, ParticleAoS<Size>>(); break;
			}

			mpi::NeighborFinder<Size>* nfp;

			switch (nf)
			{
			case 0: nfp = new mpi::SpatialHashing<Size>(); break;
			default: nfp = new mpi::MortonSorting<Size>(initp->GetDomain().first, initp->GetDomain().second); break;
			}

			mpi::SPHSimulation<Size> sim(nfp);

			if (nf == 0)
				sim.SetName("MPI AoS Hashing");
			else
				sim.SetName("MPI AoS Morton");

			sim.SetRank(mpi_rank, mpi_comm, mpi_size);
			if (mpi_rank == 0)
			{
				sim.InitializeFluid(initp);
				
				if (file)
				{
					exporter.SetBaseName(filename);
					exporter.SetFormat(fileformat);
					sim.AddObserver(&exporter);
				}
			}

			sim.Start();

			delete nfp;
			delete initp;

			// Finalize MPI.
			MPI_Finalize();
		}

		else
		{

			if (layout == 0)
			{
				using Particles = ParticleAoS<Size>;
				FileExporter<Size, Particles> exporter;

				SimInitializer<Size, Particles>* initp;

				switch (init)
				{
				default: initp = new TrayInitializer<Size, Particles>(); break;
				}

				if (par == 0)
				{
					serial::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new serial::SpatialHashing<Size, Particles>(); break;
					default: nfp = new serial::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					serial::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("Serial AoS Hashing");
					else
						sim.SetName("Serial AoS Morton");

					sim.InitializeFluid(initp);

					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				else if (par == 1)
				{
					openmp::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new openmp::SpatialHashing<Size, Particles>(); break;
					default: nfp = new openmp::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					openmp::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("OMP AoS Hashing");
					else
						sim.SetName("OMP AoS Morton");

					sim.InitializeFluid(initp);

					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				delete initp;
			}

			else if (layout == 1)
			{
				using Particles = ParticleSoA<Size>;
				FileExporter<Size, Particles> exporter;

				SimInitializer<Size, Particles>* initp;

				switch (init)
				{
				default: initp = new TrayInitializer<Size, Particles>(); break;
				}

				if (par == 0)
				{
					serial::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new serial::SpatialHashing<Size, Particles>(); break;
					default: nfp = new serial::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					serial::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("Serial SoA Hashing");
					else
						sim.SetName("Serial SoA Morton");

					sim.InitializeFluid(initp);

					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				else if (par == 1)
				{
					openmp::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new openmp::SpatialHashing<Size, Particles>(); break;
					default: nfp = new openmp::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					openmp::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("OMP SoA Hashing");
					else
						sim.SetName("OMP SoA Morton");

					sim.InitializeFluid(initp);

					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				delete initp;
			}

			else if (layout == 2)
			{
				using Particles = ParticleHybrid<Size>;
				FileExporter<Size, Particles> exporter;

				SimInitializer<Size, Particles>* initp;

				switch (init)
				{
				default: initp = new TrayInitializer<Size, Particles>(); break;
				}

				if (par == 0)
				{
					serial::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new serial::SpatialHashing<Size, Particles>(); break;
					default: nfp = new serial::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					serial::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("Serial Hybrid Hashing");
					else
						sim.SetName("Serial Hybrid Morton");

					sim.InitializeFluid(initp);

					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				else if (par == 1)
				{
					openmp::NeighborFinder<Size, Particles>* nfp;

					switch (nf)
					{
					case 0: nfp = new openmp::SpatialHashing<Size, Particles>(); break;
					default: nfp = new openmp::MortonSorting<Size, Particles>(initp->GetDomain().first, initp->GetDomain().second); break;
					}

					openmp::SPHSimulation<Size, Particles> sim(nfp);

					if (nf == 0)
						sim.SetName("OMP Hybrid Hashing");
					else
						sim.SetName("OMP Hybrid Morton");

					sim.InitializeFluid(initp);

					if (file)
					{
						exporter.SetBaseName(filename);
						exporter.SetFormat(fileformat);
						sim.AddObserver(&exporter);
					}

					sim.Start();

					delete nfp;
				}

				delete initp;
			}
		}
	}
}