# AMSC_SPH


## Dependencies

	sudo apt install wayland-protocols libwayland-dev libxkbcommon-dev pkg-config libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl1-mesa-dev mesa-common-dev openmpi-bin libopenmpi-dev

## Build

	cmake -S . --preset <windows | linux>-<debug | release | reldeb>
	cd build/<debug | release | reldeb>
	make -j

	./AMSC_SPH [options]

Use the following cmake flag to exclude UI files from the compilation

	-DDISABLE_UI=ON

## Running Examples

It's possible to execute the already available simulations using the following command line parameters
- `-example`
- `-d <2 | 3>`: simulation dimension, either 2D or 3D
- `-par <serial | openmp | mpi>`: parallelization method
        - `-init <box | tray | TODO>`: initial state of the particles, see [Defining Custom Simulation](#defining-custom-simulation) for further details
- `-nf <hash | morton>`: strategy to find particle neighbors, either spatial hashing or morton sorting
- `-layout <AoS | SoA | hybrid>`: memory layout for collections of particles, see [Particles Memory Layout](#particles-memory-layout) for more info
- `[-file <filename> <xyz | vtu | bvtu>]`: if present, writes particle state to a file with format ".xyz", ".vtu" or ".vtu" with binary encoding, and name `filename-{frame number}`
- `[-ui]`: if present, use realtime visualization (only for 2D)
- `[--omp-threads <num>]`:
	using the number of physical cores might yield a faster execution,
	by default OpenMP uses the number of logical / SMT / Hyper-Threading cores,
	which is worse because SPH can be memory-bandwidth bound, not compute-bound

## Defining Custom Simulation

Each simulation is defined by the following elements
- dimension: either 2D or 3D
- particles memory layout
- parallelization method
- neighbor finding method
- simulation parameters
- initial state of the particles
- simulation progress observer

### Particles Memory Layout

The code defines three different memory layouts for a collection of particles
- ParticleAoS: each particle data is contiguous in memory
- ParticleSoA: defines separate vectors for all the atomic informations that describe a particle, 
			   this means that there's a vector for the x component of the position of each particle,
			   followed by a vector for all y components, and so on
- ParticleHybrid: intermediate design where vector (in the mathematical sense) components are not separated in distinct std::vector

> Note: MPI parallelization supports only ParticleAoS, support for other memory layouts would have required a specific implementation for each of them

### Parallelization Method

There are three available parallelizations: serial, OpenMP and MPI.

Each of them has a corresponding namespace (`serial`, `openmp`, `mpi`) and
a corresponding implementation of the class `SPHSimulation`, extending and overriding
methods from the base class `base::SPHSimulation`.

Each parallelization method also defines custom implementations of the neighbor finding strategies.
See [Neighbor Finding Method](#neighbor-finding-method)

> Note: a simulation with a certain parallelization cannot be combined with a neighbor finding strategy
with a different parallelization, as they make different assumptions on how they are used.

> Note: in the following sections we will use the `par` namespace to refer to either
`serial`, `openmp` or `mpi`, and `sim` will refer to an instance of `par::SHPSimulation`

### Neighbor Finding Method

There are currently two different strategies implemented: `par::SpatialHashing` and `par::MortonSorting`; 
the latter yielding better performance. Both are templated to support different memory layouts.

An instance of the desired method has to be passed to the simulation via constructor argument.

> Note: the memory layout used in the neighbor finding method has to match that of the simulation

### Simulation Parameters

The available simulation parameters relevant to the system evolution are those present in the struct `base::SPHParameters`
and they include timestep or physical particle properties.
The simulation parameters can be specified via the function `sim.SetParams(params)`.

They can be updated during the execution of the simulation, between frames. See [Simulation Progress Observer](#simulation-progress-observer) section.

### Initial State of the Particles

The currently available initializers are
- 2D `BoxInitializer`: creates a rectangle of solid particles with fluid particles inside, it also adds some solid obstacles
- 3D `TrayInitializer`: it's like a 3D box sliced horizontally, keeping the lower half.
- TODO: add the other initializers

Both classes extend the interface `SimInitializer`, and require the memory layout has template argument.
A custom initializer can be defined by extending `SimInitializer` and implementing the `Init` method,
and making use of the utility functions to add walls, obstacles, sources and sinks of particles.

At the start of a simulation the particles have to be initialized by calling the function `sim.InitializeFluid(&initializer)`.

> Note: the memory layout used in the initializer has to match that of the simulation

### Simulation Progress Observer

An observer for the simulation is a class that inherits from the base class `Observer`.
It's possible to attach observers to the simulation using the method `sim.AddObserver(&observer)`.
The simulation will notify the attached observers of the start and end of each frame, calling the methods
`observer->OnStartFrame()` and `observer->OnEndFrame()`

During either of them it's safe to update or fetch the simulation parameters and particle data.

The currently implemented observers are:
- FileExporter: this allows writing the particles state to ".xyz" or ".vtu" files; given the small timestep required
	by some simulations, it allows actually saving to file once every fixed amount of frames.
- ImGuiViewer: realtime visualization of 2D simulations, also allows altering simulation parameters
	and applying pressure to the fluid by clicking it; efficient up to simulations in the order of 100'000 particles

> Note: ImGuiViewer requires compiling without the -DDISABLE_UI=ON cmake flag

### Putting Everything Together

Below there is a full example of initialization

	constexpr size_t Size = 2;
	using Particles = ParticleAoS<Size>;

	openmp::SpatialHashing<Size, Particles> nf;

	openmp::SPHSimulation<Size, Particles> sph(&nf);
	sph.SetName("AoS-OMP-Hashing");

	BoxInitializer<Size, Particles> initializer;
	sph.InitializeFluid(&initializer);

	ImGuiViewer<Particles> imguiViewer;
	sph.AddObserver(&imguiViewer);
	imguiViewer.Start();
	sph.Start();