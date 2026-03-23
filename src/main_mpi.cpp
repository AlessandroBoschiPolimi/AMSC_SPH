#include <mpi.h>
#include <memory>
#include "Utility.hpp"
#include "MPI/SPHSimulation.hpp"

#include "Visualizers/FileExporter.hpp"
#include "Visualizers/ImGuiViewer.hpp"

#include "MPI/Neighbors/SpatialHashing.hpp"
#include "MPI/Neighbors/MortonSorting.hpp"

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"
int main( int argc , char **argv)
{
MPI_Init(&argc , &argv );
MPI_Comm mpi_comm = MPI_COMM_WORLD;
//Get the number of processes .
int mpi_size;
MPI_Comm_size(mpi_comm, &mpi_size );
//Get the rank of current process .
int mpi_rank;
MPI_Comm_rank(mpi_comm, &mpi_rank);
{
	constexpr size_t Size = 3;
	using Particles = ParticleAoS<Size>;

	//mpi::SpatialHashing<Size, Particles> nf;
	mpi::MortonSorting<Size, Particles> nf({ 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f });
	mpi::SPHSimulation<Size, Particles> sph(&nf);
	sph.SetRank(mpi_rank, mpi_comm, mpi_size);
	if (mpi_rank == 0)
	{
		FileExporter<Size, Particles> imguiViewer;
		imguiViewer.SetFrequency(20);
		TrayInitializer<Size, Particles> initializer;
		sph.InitializeFluid(&initializer);

		sph.AddObserver(&imguiViewer);
		//imguiViewer.Start();
		sph.Start();
	}
	else
		sph.Start();
}
// Finalize MPI.
MPI_Finalize();
}
