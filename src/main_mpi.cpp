#include <mpi.h>
#include <memory>
#include "Utility.hpp"
#include "MPI/SPHSimulation.hpp"

#include "Visualizers/FileExporter.hpp"

#ifndef DISABLE_UI
#include "Visualizers/ImGuiViewer.hpp"
#endif

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
	constexpr size_t Size = 2;

	//mpi::SpatialHashing<Size, Particles> nf;
	//mpi::MortonSorting<Size> nf({ 0.0f, 0.0f}, { 1.0f, 1.0f});
	mpi::SpatialHashing<Size> nf;
	mpi::SPHSimulation<Size> sph(&nf);
	sph.SetRank(mpi_rank, mpi_comm, mpi_size);
	if (mpi_rank == 0)
	{
		using Particles = ParticleAoS<Size>;
#ifndef DISABLE_UI
		ImGuiViewer<Particles> imguiViewer;
#endif
		BoxInitializer<Size, Particles> initializer;
		sph.InitializeFluid(&initializer);

#ifndef DISABLE_UI
		sph.AddObserver(&imguiViewer);
		imguiViewer.Start();
#endif
		sph.Start();
	}
	else
		sph.Start();
}
// Finalize MPI.
MPI_Finalize();
}
