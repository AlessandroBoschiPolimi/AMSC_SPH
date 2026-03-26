#pragma once
#include "Particle.hpp"
#include "Base/SPHSimulation.hpp"

#include "OpenMP/SPHSimulation.hpp"
#include "Serial/SPHSimulation.hpp"

#include "Visualizers/FileExporter.hpp"
#include "Visualizers/ImGuiViewer.hpp"

#include "Initializers/TrayInitializer.hpp"
#include "Initializers/BoxInitializer.hpp"

#include "OpenMP/Neighbors/SpatialHashing.hpp"
#include "OpenMP/Neighbors/MortonSorting.hpp"

#include "Serial/Neighbors/SpatialHashing.hpp"
#include "Serial/Neighbors/MortonSorting.hpp"


void Test_2D_OpenMP_Morton();
void Test_2D_OpenMP_Hashing();
void Test_2D_Serial_Morton();
void Test_2D_Serial_Hashing();
void Test_2D_Pipe_OpenMP_Morton();
void Test_Correctness_Generate();
void Test_Correctness_Check();

/// Dont call the following functions directly.

template <size_t D, ParticleSet<D> Particles>
serial::MortonSorting<D, Particles> GenerateTest_serial_MortonSorting(const SimInitializer<D, Particles>* i)
{
	if constexpr (D == 2)
		return serial::MortonSorting<D, Particles>(i->GetDomain().first, i->GetDomain().second);
	else 
		return serial::MortonSorting<D, Particles>(i->GetDomain().first, i->GetDomain().second);
}
template <size_t D, ParticleSet<D> Particles>
openmp::MortonSorting<D, Particles> GenerateTest_openmp_MortonSorting(const SimInitializer<D, Particles>* i)
{
	if constexpr (D == 2)
		return openmp::MortonSorting<D, Particles>(i->GetDomain().first, i->GetDomain().second);
	else
		return openmp::MortonSorting<D, Particles>(i->GetDomain().first, i->GetDomain().second);
}
template <size_t D, ParticleSet<D> Particles>
serial::SpatialHashing<D, Particles> GenerateTest_serial_SpatialHashing(const SimInitializer<D, Particles>* i)
{ return serial::SpatialHashing<D, Particles>{}; }
template <size_t D, ParticleSet<D> Particles>
openmp::SpatialHashing<D, Particles> GenerateTest_openmp_SpatialHashing(const SimInitializer<D, Particles>* i)
{ return openmp::SpatialHashing<D, Particles>{}; }

template <size_t D, ParticleSet<D> Particles>
consumer<ImGuiViewer<Particles>*> GenerateTest_View()
{ return [](ImGuiViewer<Particles>* v) { v->Start(); }; }
template <size_t D, ParticleSet<D> Particles>
consumer<FileExporter<D, Particles>*> GenerateTest_View(const std::string& path, const ExportFormat& format)
{
	return [&](FileExporter<D, Particles>* v) {
			v->SetBaseName(path);
			v->SetFormat(format);
		};
}

/// In the followin macros
/// - s: size, int
/// - ns: namespace, serial/openmp/...
/// - nf: neighbor finder, SpatialHashing/MortonSorting
/// - particles: class that satisfies ParticleSet<s>, ParticleAoS/ParticleSoA/...
/// - view: Observer, ImGuiViewer/FileExporter
/// - init: SimInitializer, BoxInitializer/...
///
/// Overall those macros are just to make testing less tedious, reducing verbosity of calling "Test" function, eg:
///		GENERATE_TEST_CASE_VIEW(2, openmp, MortonSorting, ParticleAoS, BoxInitializer, "2D_OMP_AoS_Morton", base::SPHParams{});
/// becomes
///		Test
///			<2, openmp::MortonSorting<2, ParticleAoS<2>>, openmp::SPHSimulation<2, ParticleAoS<2>>, ImGuiViewer<ParticleAoS<2>>, BoxInitializer<2, ParticleAoS<2>>, ParticleAoS<2>>
///			("2D_OMP_AoS_Morton", base::SPHParams{}, []() -> openmp::MortonSorting<2, ParticleAoS<2>> { return GenerateTest_openmp_MortonSorting<2, ParticleAoS<2>>(); }, GenerateTest_View<2, ParticleAoS<2>>());
/// 
/// they shouldn't be used for anything else

#define GENERATE_TEST_TEMPLATE(s, ns, nf, particles, view, init) \
	s, ns::nf<s, particles<s>>, ns::SPHSimulation<s, particles<s>>, view, init<s, particles<s>>, particles<s>
#define GENERATE_FILE_EXPORTER_TYPE(s, particles) FileExporter<s, particles<s>>
#define GENERATE_TEST_NF(s, ns, nf, particles) \
	[](const SimInitializer<s, particles<s>>* i) -> ns::nf<s, particles<s>> \
	{ return GenerateTest_##ns##_##nf<s, particles<s>>(i); }

/// use to launch a visualization 
#define GENERATE_TEST_CASE_VIEW(s, ns, nf, particles, init, name, params) \
	Test<GENERATE_TEST_TEMPLATE(s, ns, nf, particles, ImGuiViewer<particles<s>>, init)>( \
		name, params, \
		GENERATE_TEST_NF(s, ns, nf, particles), \
		GenerateTest_View<s, particles<s>>() \
	)
/// use to launch an execution that produces a file for comparison
#define GENERATE_TEST_CASE_COMPARE(s, ns, nf, particles, output, init, name, params) \
	Test<GENERATE_TEST_TEMPLATE(s, ns, nf, particles, GENERATE_FILE_EXPORTER_TYPE(s, particles), init)>( \
		name, params, \
		GENERATE_TEST_NF(s, ns, nf, particles), \
		GenerateTest_View<s, particles<s>>(output, LANDMINE) \
	)
/// use to launch an execution that produces a VTU file
#define GENERATE_TEST_CASE_PARAVIEW(s, ns, nf, particles, output, init, name, params) \
	Test<GENERATE_TEST_TEMPLATE(s, ns, nf, particles, GENERATE_FILE_EXPORTER_TYPE(s, particles), init)>( \
		name, params, \
		GENERATE_TEST_NF(s, ns, nf, particles), \
		GenerateTest_View<s, particles<s>>(output, VTU_01) \
	)


template <size_t D, typename NF, typename SIM, typename VIEW, typename INIT, ParticleSet<D> Particles>
void Test(const std::string& name, const base::SPHParams& params, const std::function<NF(const INIT* init)>& init_neigh, const consumer<VIEW*>& init_view)
{
	INIT init;
	NF nf = init_neigh(&init);

	SIM sph(&nf);
	sph.SetName(name);
	sph.SetParams(params);
	sph.InitializeFluid(&init);

	VIEW viewer;
	sph.AddObserver(&viewer);
	init_view(&viewer);

	sph.Start();
}
