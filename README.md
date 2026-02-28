# AMSC_SPH


## Dependencies

	sudo apt install wayland-protocols libwayland-dev libxkbcommon-dev pkg-config libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl1-mesa-dev mesa-common-dev openmpi-bin libopenmpi-dev

## Build

	cmake -S . --preset <windows | linux>-<debug | release | reldeb>

	cd build/<debug | release | reldeb>
	make -j

	./AMSC_SPH [options]

Use the following cmake flag to exclude cuda files from the compilation

	-DENABLE_CUDA=OFF


### options
- `--omp-threads <num>`: 
	using the number of physical cores might yield a faster execution,
	by default OpenMP uses the number of logical / SMT / Hyper-Threading cores,
	which is worse because SPH can be memory-bandwidth bound, not compute-bound