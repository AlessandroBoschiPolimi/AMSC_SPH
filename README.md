# AMSC_SPH


## Dependencies

	sudo apt install wayland-protocols libwayland-dev libxkbcommon-dev pkg-config libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl1-mesa-dev mesa-common-dev

## Build

	cmake -S . --preset <windows | linux>-<debug | release>

	cd build/<debug | release>
	make -j

	./AMSC_SPH [options]

### options
- `--omp-threads <num>`: 
	using the number of physical cores is suggested,
	by default OpenMP uses the number of logical / SMT / Hyper-Threading cores,
	which is worse because SPH is memory-bandwidth bound, not compute-bound
