#pragma once
#include "Utility.hpp"

enum ParticleType : u8
{
	FLUID,
	SOLID
};

template <size_t D>
struct Particle
{
	using vec_t = coord<float, D>;
	using idx_t = u32; // alias to specify that the int represents the index of a particles

	vec_t Position;
	vec_t Velocity;
	vec_t A_grav;
	vec_t A_visc;
	vec_t A_press;
	float Mass;
	float Density;
	float Pressure;
	ParticleType Type;
	float BoundaryPsi;
};


/*
AoS   : [{x, y, z, vx, vy, vz, ...}]             Array of particles
Hybrid: {[{x, y, z}], [{vx, vy, vz}], ...}       One array per coord
SoA   : {[x], [y], [z], [vx], [vy], [vz], ...}   One array per component
*/

template <typename P, size_t D>
concept ParticleSet = std::movable<P> && std::default_initializable<P> && 
requires(const P cp, P p, P other, size_t i, const Particle<D>& particle, Particle<D>&& particle_move, std::vector<Particle<D>>& out)
{
	{ cp.Size()				} -> std::convertible_to<size_t>;
	{ cp.Empty()			} -> std::convertible_to<bool>;
	{ p.Resize(size_t{})    };
	{ p.Clear()				};
	{ p.PushBack(particle)	};

	{ cp.GetParticle(i)					} -> std::same_as<Particle<D>>;
	{ cp.GetParticles(out)				};
	{ p.SetParticle(i, particle)		};
	{ p.SetParticle(i, particle_move)	};
	{ p.SetParticle(i, other, size_t{}) };

	{ cp.Position(i)				} -> std::convertible_to<typename Particle<D>::vec_t>;
	{ cp.PositionX(i)				} -> std::convertible_to<float>;
	{ cp.PositionY(i)				} -> std::convertible_to<float>;
	{ cp.PositionZ(i)				} -> std::convertible_to<float>;
	
	{ cp.Velocity(i)				} -> std::convertible_to<typename Particle<D>::vec_t>;
	{ cp.VelocityX(i)				} -> std::convertible_to<float>;
	{ cp.VelocityY(i)				} -> std::convertible_to<float>;
	{ cp.VelocityZ(i)				} -> std::convertible_to<float>;

	{ cp.AccelerationGravity(i)		} -> std::convertible_to<typename Particle<D>::vec_t>;
	{ cp.AccelerationGravityX(i)	} -> std::convertible_to<float>;
	{ cp.AccelerationGravityY(i)	} -> std::convertible_to<float>;
	{ cp.AccelerationGravityZ(i)	} -> std::convertible_to<float>;

	{ cp.AccelerationViscosity(i)	} -> std::convertible_to<typename Particle<D>::vec_t>;
	{ cp.AccelerationViscosityX(i)	} -> std::convertible_to<float>;
	{ cp.AccelerationViscosityY(i)	} -> std::convertible_to<float>;
	{ cp.AccelerationViscosityZ(i)	} -> std::convertible_to<float>;

	{ cp.AccelerationPressure(i)	} -> std::convertible_to<typename Particle<D>::vec_t>;
	{ cp.AccelerationPressureX(i)	} -> std::convertible_to<float>;
	{ cp.AccelerationPressureY(i)	} -> std::convertible_to<float>;
	{ cp.AccelerationPressureZ(i)	} -> std::convertible_to<float>;

	{ cp.Mass(i) 					} -> std::convertible_to<float>;
	{ cp.Density(i) 				} -> std::convertible_to<float>;
	{ cp.Pressure(i) 				} -> std::convertible_to<float>;
	{ cp.Type(i) 					} -> std::convertible_to<ParticleType>;
	{ cp.BoundaryPsi(i) 			} -> std::convertible_to<float>;

	{ p.SetPosition(i, typename Particle<D>::vec_t{}) };
	{ p.SetPositionX(i, float{}) };
	{ p.SetPositionY(i, float{}) };
	{ p.SetPositionZ(i, float{}) };

	{ p.SetVelocity(i, typename Particle<D>::vec_t{}) };
	{ p.SetVelocityX(i, float{}) };
	{ p.SetVelocityY(i, float{}) };
	{ p.SetVelocityZ(i, float{}) };
	
	{ p.SetAccelerationGravity(i, typename Particle<D>::vec_t{}) };
	{ p.SetAccelerationGravityX(i, float{}) };
	{ p.SetAccelerationGravityY(i, float{}) };
	{ p.SetAccelerationGravityZ(i, float{}) };

	{ p.SetAccelerationViscosity(i, typename Particle<D>::vec_t{}) };
	{ p.SetAccelerationViscosityX(i, float{}) };
	{ p.SetAccelerationViscosityY(i, float{}) };
	{ p.SetAccelerationViscosityZ(i, float{}) };

	{ p.SetAccelerationPressure(i, typename Particle<D>::vec_t{}) };
	{ p.SetAccelerationPressureX(i, float{}) };
	{ p.SetAccelerationPressureY(i, float{}) };
	{ p.SetAccelerationPressureZ(i, float{}) };

	{ p.SetMass(i, float{})			};
	{ p.SetDensity(i, float{})		};
	{ p.SetPressure(i, float{})		};
	{ p.SetType(i, ParticleType{})	};
	{ p.SetBoundaryPsi(i, float{})	};
};


/// Calling Z component functions when D == 2 is undefined behaviour
template <size_t D>
struct ParticleHybrid
{
	using vec_t = Particle<D>::vec_t;
	using idx_t = Particle<D>::idx_t;

	std::vector<vec_t> Positions;
	std::vector<vec_t> Velocitys;
	std::vector<vec_t> A_gravs;
	std::vector<vec_t> A_viscs;
	std::vector<vec_t> A_presss;
	std::vector<float> Masss;
	std::vector<float> Densitys;
	std::vector<float> Pressures;
	std::vector<ParticleType> Types;
	std::vector<float> BoundaryPsis;

	ParticleHybrid() = default;
	ParticleHybrid(ParticleHybrid&& other)
		: Positions(std::move(other.Positions)), Velocitys(std::move(other.Velocitys)), A_gravs(std::move(other.A_gravs))
		, A_viscs(std::move(other.A_viscs)), A_presss(std::move(other.A_presss)), Masss(std::move(other.Masss))
		, Densitys(std::move(other.Densitys)), Pressures(std::move(other.Pressures)), Types(std::move(other.Types))
		, BoundaryPsis(std::move(other.BoundaryPsis))
	{}

	ParticleHybrid& operator=(const ParticleHybrid& other) {
		Positions = other.Positions;
		Velocitys = other.Velocitys;
		A_gravs = other.A_gravs;
		A_viscs = other.A_viscs;
		A_presss = other.A_presss;
		Masss = other.Masss;
		Densitys = other.Densitys;
		Pressures = other.Pressures;
		Types = other.Types;
		BoundaryPsis = other.BoundaryPsis;
		return *this;
	}
	ParticleHybrid& operator=(ParticleHybrid&& other) noexcept {
		Positions = std::move(other.Positions);
		Velocitys = std::move(other.Velocitys);
		A_gravs = std::move(other.A_gravs);
		A_viscs = std::move(other.A_viscs);
		A_presss = std::move(other.A_presss);
		Masss = std::move(other.Masss);
		Densitys = std::move(other.Densitys);
		Pressures = std::move(other.Pressures);
		Types = std::move(other.Types);
		BoundaryPsis = std::move(other.BoundaryPsis);
		return *this;
	}


	size_t Size() const { return Positions.size(); }
	bool Empty() const { return Positions.empty(); }
	void Resize(size_t size) {
		Positions.resize(size);
		Velocitys.resize(size);
		A_gravs.resize(size);
		A_viscs.resize(size);
		A_presss.resize(size);
		Masss.resize(size);
		Densitys.resize(size);
		Pressures.resize(size);
		Types.resize(size);
		BoundaryPsis.resize(size);
	}
	void Clear() {
		Positions.clear();
		Velocitys.clear();
		A_gravs.clear();
		A_viscs.clear();
		A_presss.clear();
		Masss.clear();
		Densitys.clear();
		Pressures.clear();
		Types.clear();
		BoundaryPsis.clear();
	}
	void PushBack(const Particle<D>& particle) {
		Positions.push_back(particle.Position);
		Velocitys.push_back(particle.Velocity);
		A_gravs.push_back(particle.A_grav);
		A_viscs.push_back(particle.A_visc);
		A_presss.push_back(particle.A_press);
		Masss.push_back(particle.Mass);
		Densitys.push_back(particle.Density);
		Pressures.push_back(particle.Pressure);
		Types.push_back(particle.Type);
		BoundaryPsis.push_back(particle.BoundaryPsi);
	}

	Particle<D> GetParticle(size_t i) const {
		Particle<D> p;
		p.Position    = Positions[i];
		p.Velocity    = Velocitys[i];
		p.A_grav      = A_gravs[i];
		p.A_visc      = A_viscs[i];
		p.A_press     = A_presss[i];
		p.Mass        = Masss[i];
		p.Density     = Densitys[i];
		p.Pressure    = Pressures[i];
		p.Type        = Types[i];
		p.BoundaryPsi = BoundaryPsis[i];
		return p;
	}
	void GetParticles(std::vector<Particle<D>>& out) const {
		out.clear();
		out.resize(Positions.size());
		
		for (size_t i = 0; i < Positions.size(); ++i)
			out[i] = GetParticle(i);
	}
	void SetParticle(size_t i, const Particle<D>& particle) {
		Positions[i] = particle.Position;
		Velocitys[i] = particle.Velocity;
		A_gravs[i] = particle.A_grav;
		A_viscs[i] = particle.A_visc;
		A_presss[i] = particle.A_press;
		Masss[i] = particle.Mass;
		Densitys[i] = particle.Density;
		Pressures[i] = particle.Pressure;
		Types[i] = particle.Type;
		BoundaryPsis[i] = particle.BoundaryPsi;
	}
	void SetParticle(size_t i, Particle<D>&& particle) {
		Positions[i] = std::move(particle.Position);
		Velocitys[i] = std::move(particle.Velocity);
		A_gravs[i] = std::move(particle.A_grav);
		A_viscs[i] = std::move(particle.A_visc);
		A_presss[i] = std::move(particle.A_press);
		Masss[i] = particle.Mass;
		Densitys[i] = particle.Density;
		Pressures[i] = particle.Pressure;
		Types[i] = particle.Type;
		BoundaryPsis[i] = particle.BoundaryPsi;
	}
	void SetParticle(size_t i, const ParticleHybrid& other, size_t src) {
		Positions[i] = other.Positions[src];
		Velocitys[i] = other.Velocitys[src];
		A_gravs[i] = other.A_gravs[src];
		A_viscs[i] = other.A_viscs[src];
		A_presss[i] = other.A_presss[src];
		Masss[i] = other.Masss[src];
		Densitys[i] = other.Densitys[src];
		Pressures[i] = other.Pressures[src];
		Types[i] = other.Types[src];
		BoundaryPsis[i] = other.BoundaryPsis[src];
	}

	vec_t Position (size_t i) const { return Positions[i]; }
	float PositionX(size_t i) const { return Positions[i].x; }
	float PositionY(size_t i) const { return Positions[i].y; }
	float PositionZ(size_t i) const { if constexpr (D > 2) return Positions[i].z; else return 0.f; }

	vec_t Velocity (size_t i) const { return Velocitys[i]; }
	float VelocityX(size_t i) const { return Velocitys[i].x; }
	float VelocityY(size_t i) const { return Velocitys[i].y; }
	float VelocityZ(size_t i) const { if constexpr (D > 2) return Velocitys[i].z; else return 0.f; }

	vec_t AccelerationGravity (size_t i) const { return A_gravs[i]; }
	float AccelerationGravityX(size_t i) const { return A_gravs[i].x; }
	float AccelerationGravityY(size_t i) const { return A_gravs[i].y; }
	float AccelerationGravityZ(size_t i) const { if constexpr (D > 2) return A_gravs[i].z; else return 0.f; }

	vec_t AccelerationViscosity (size_t i) const { return A_viscs[i]; }
	float AccelerationViscosityX(size_t i) const { return A_viscs[i].x; }
	float AccelerationViscosityY(size_t i) const { return A_viscs[i].y; }
	float AccelerationViscosityZ(size_t i) const { if constexpr (D > 2) return A_viscs[i].z; else return 0.f; }

	vec_t AccelerationPressure(size_t i)  const { return A_presss[i]; }
	float AccelerationPressureX(size_t i) const { return A_presss[i].x; }
	float AccelerationPressureY(size_t i) const { return A_presss[i].y; }
	float AccelerationPressureZ(size_t i) const { if constexpr (D > 2) return A_presss[i].z; else return 0.f; }

	float Mass       (size_t i) const { return Masss[i]; }
	float Density    (size_t i) const { return Densitys[i]; }
	float Pressure   (size_t i) const { return Pressures[i]; }
	ParticleType Type(size_t i) const { return Types[i]; }
	float BoundaryPsi(size_t i) const { return BoundaryPsis[i]; }

	void SetPosition(size_t i, const vec_t& pos) { Positions[i] = pos; }
	void SetPositionX(size_t i, float v) { Positions[i].x = v; }
	void SetPositionY(size_t i, float v) { Positions[i].y = v; }
	void SetPositionZ(size_t i, float v) { if constexpr (D > 2) Positions[i].z = v; }

	void SetVelocity(size_t i, const vec_t& vel) { Velocitys[i] = vel; }
	void SetVelocityX(size_t i, float v) { Velocitys[i].x = v; }
	void SetVelocityY(size_t i, float v) { Velocitys[i].y = v; }
	void SetVelocityZ(size_t i, float v) { if constexpr (D > 2) Velocitys[i].z = v; }

	void SetAccelerationGravity(size_t i, const vec_t& acc) { A_gravs[i] = acc; }
	void SetAccelerationGravityX(size_t i, float v) { A_gravs[i].x = v; }
	void SetAccelerationGravityY(size_t i, float v) { A_gravs[i].y = v; }
	void SetAccelerationGravityZ(size_t i, float v) { if constexpr (D > 2) A_gravs[i].z = v; }

	void SetAccelerationViscosity(size_t i, const vec_t& acc) { A_viscs[i] = acc; }
	void SetAccelerationViscosityX(size_t i, float v) { A_viscs[i].x = v; }
	void SetAccelerationViscosityY(size_t i, float v) { A_viscs[i].y = v; }
	void SetAccelerationViscosityZ(size_t i, float v) { if constexpr (D > 2) A_viscs[i].z = v; }

	void SetAccelerationPressure(size_t i, const vec_t& acc) { A_presss[i] = acc; }
	void SetAccelerationPressureX(size_t i, float v) { A_presss[i].x = v; }
	void SetAccelerationPressureY(size_t i, float v) { A_presss[i].y = v; }
	void SetAccelerationPressureZ(size_t i, float v) { if constexpr (D > 2) A_presss[i].z = v; }

	void SetMass(size_t i, float m) { Masss[i] = m; }
	void SetDensity(size_t i, float d) { Densitys[i] = d; }
	void SetPressure(size_t i, float p) { Pressures[i] = p; }
	void SetType(size_t i, ParticleType t) { Types[i] = t; }
	void SetBoundaryPsi(size_t i, float p) { BoundaryPsis[i] = p; }
};


/// Calling Z component functions when D == 2 is undefined behaviour
template <size_t D>
struct ParticleSoA
{
	using vec_t = typename Particle<D>::vec_t;
	using idx_t = typename Particle<D>::idx_t;

	std::vector<float> Xs, Ys, Zs;
	std::vector<float> VXs, VYs, VZs;
	std::vector<float> AX_gravs, AY_gravs, AZ_gravs;
	std::vector<float> AX_viscs, AY_viscs, AZ_viscs;
	std::vector<float> AX_presss, AY_presss, AZ_presss;
	std::vector<float> Masss;
	std::vector<float> Densitys;
	std::vector<float> Pressures;
	std::vector<ParticleType> Types;
	std::vector<float> BoundaryPsis;


	ParticleSoA() = default;
	ParticleSoA(ParticleSoA&& other)
		: Xs(std::move(other.Xs)), Ys(std::move(other.Ys)), Zs(std::move(other.Zs))
		, VXs(std::move(other.VXs)), VYs(std::move(other.VYs)), VZs(std::move(other.VZs))
		, AX_gravs(std::move(other.AX_gravs)), AY_gravs(std::move(other.AY_gravs)), AZ_gravs(std::move(other.AZ_gravs))
		, AX_viscs(std::move(other.AX_viscs)), AY_viscs(std::move(other.AY_viscs)), AZ_viscs(std::move(other.AZ_viscs))
		, AX_presss(std::move(other.AX_presss)), AY_presss(std::move(other.AY_presss)), AZ_presss(std::move(other.AZ_presss))
		, Masss(std::move(other.Masss)), Densitys(std::move(other.Densitys)), Pressures(std::move(other.Pressures))
		, Types(std::move(other.Types)), BoundaryPsis(std::move(other.BoundaryPsis))
	{ }

	ParticleSoA& operator=(const ParticleSoA& other) {
		Xs = other.Xs; Ys = other.Ys; if constexpr (D > 2) Zs = other.Zs;
		VXs = other.VXs; VYs = other.VYs; if constexpr (D > 2) VZs = other.VZs;
		AX_gravs = other.AX_gravs; AY_gravs = other.AY_gravs; if constexpr (D > 2) AZ_gravs = other.AZ_gravs;
		AX_viscs = other.AX_viscs; AY_viscs = other.AY_viscs; if constexpr (D > 2) AZ_viscs = other.AZ_viscs;
		AX_presss = other.AX_presss; AY_presss = other.AY_presss; if constexpr (D > 2) AZ_presss = other.AZ_presss;
		Masss = other.Masss;
		Densitys = other.Densitys;
		Pressures = other.Pressures;
		Types = other.Types;
		BoundaryPsis = other.BoundaryPsis;
		return *this;
	}
	ParticleSoA& operator=(ParticleSoA&& other) noexcept {
		Xs = std::move(other.Xs); Ys = std::move(other.Ys); if constexpr (D > 2) Zs = std::move(other.Zs);
		VXs = std::move(other.VXs); VYs = std::move(other.VYs); if constexpr (D > 2) VZs = std::move(other.VZs);
		AX_gravs = std::move(other.AX_gravs); AY_gravs = std::move(other.AY_gravs); if constexpr (D > 2) AZ_gravs = std::move(other.AZ_gravs);
		AX_viscs = std::move(other.AX_viscs); AY_viscs = std::move(other.AY_viscs); if constexpr (D > 2) AZ_viscs = std::move(other.AZ_viscs);
		AX_presss = std::move(other.AX_presss); AY_presss = std::move(other.AY_presss); if constexpr (D > 2) AZ_presss = std::move(other.AZ_presss);
		Masss = std::move(other.Masss);
		Densitys = std::move(other.Densitys);
		Pressures = std::move(other.Pressures);
		Types = std::move(other.Types);
		BoundaryPsis = std::move(other.BoundaryPsis);
		return *this;
	}


	// -------------------------------------------------
	// Size / Memory
	// -------------------------------------------------

	size_t Size() const { return Xs.size(); }
	bool Empty() const { return Xs.empty(); }
	void Resize(size_t n)
	{
		Xs.resize(n); Ys.resize(n); if constexpr (D > 2) Zs.resize(n);
		VXs.resize(n); VYs.resize(n); if constexpr (D > 2) VZs.resize(n);
		AX_gravs.resize(n); AY_gravs.resize(n); if constexpr (D > 2) AZ_gravs.resize(n);
		AX_viscs.resize(n); AY_viscs.resize(n); if constexpr (D > 2) AZ_viscs.resize(n);
		AX_presss.resize(n); AY_presss.resize(n); if constexpr (D > 2) AZ_presss.resize(n);
		Masss.resize(n);
		Densitys.resize(n);
		Pressures.resize(n);
		Types.resize(n);
		BoundaryPsis.resize(n);
	}
	void Clear()
	{
		Xs.clear(); Ys.clear(); Zs.clear();
		VXs.clear(); VYs.clear(); VZs.clear();
		AX_gravs.clear(); AY_gravs.clear(); AZ_gravs.clear();
		AX_viscs.clear(); AY_viscs.clear(); AZ_viscs.clear();
		AX_presss.clear(); AY_presss.clear(); AZ_presss.clear();
		Masss.clear();
		Densitys.clear();
		Pressures.clear();
		Types.clear();
		BoundaryPsis.clear();
	}

	// -------------------------------------------------
	// Conversion
	// -------------------------------------------------

	void PushBack(const Particle<D>& p)
	{
		Xs.push_back(p.Position.x);
		Ys.push_back(p.Position.y);
		if constexpr (D > 2) Zs.push_back(p.Position.z);

		VXs.push_back(p.Velocity.x);
		VYs.push_back(p.Velocity.y);
		if constexpr (D > 2) VZs.push_back(p.Velocity.z);

		AX_gravs.push_back(p.A_grav.x);
		AY_gravs.push_back(p.A_grav.y);
		if constexpr (D > 2) AZ_gravs.push_back(p.A_grav.z);

		AX_viscs.push_back(p.A_visc.x);
		AY_viscs.push_back(p.A_visc.y);
		if constexpr (D > 2) AZ_viscs.push_back(p.A_visc.z);

		AX_presss.push_back(p.A_press.x);
		AY_presss.push_back(p.A_press.y);
		if constexpr (D > 2) AZ_presss.push_back(p.A_press.z);

		Masss.push_back(p.Mass);
		Densitys.push_back(p.Density);
		Pressures.push_back(p.Pressure);
		Types.push_back(p.Type);
		BoundaryPsis.push_back(p.BoundaryPsi);
	}

	Particle<D> GetParticle(size_t i) const
	{
		Particle<D> p;

		p.Position.x = Xs[i];
		p.Position.y = Ys[i];
		if constexpr (D > 2) p.Position.z = Zs[i];

		p.Velocity.x = VXs[i];
		p.Velocity.y = VYs[i];
		if constexpr (D > 2) p.Velocity.z = VZs[i];

		p.A_grav.x = AX_gravs[i];
		p.A_grav.y = AY_gravs[i];
		if constexpr (D > 2) p.A_grav.z = AZ_gravs[i];

		p.A_visc.x = AX_viscs[i];
		p.A_visc.y = AY_viscs[i];
		if constexpr (D > 2) p.A_visc.z = AZ_viscs[i];

		p.A_press.x = AX_presss[i];
		p.A_press.y = AY_presss[i];
		if constexpr (D > 2) p.A_press.z = AZ_presss[i];

		p.Mass = Masss[i];
		p.Density = Densitys[i];
		p.Pressure = Pressures[i];
		p.Type = Types[i];
		p.BoundaryPsi = BoundaryPsis[i];

		return p;
	}
	void GetParticles(std::vector<Particle<D>>& out) const {
		out.clear();
		out.resize(Xs.size());

		for (size_t i = 0; i < Xs.size(); ++i)
			out[i] = GetParticle(i); // TODO: improve, populating all positions first...
	}

	void SetParticle(size_t i, const Particle<D>& p)
	{
		Xs[i] = p.Position.x;
		Ys[i] = p.Position.y;
		if constexpr (D > 2) Zs[i] = p.Position.z;

		VXs[i] = p.Velocity.x;
		VYs[i] = p.Velocity.y;
		if constexpr (D > 2) VZs[i] = p.Velocity.z;

		AX_gravs[i] = p.A_grav.x;
		AY_gravs[i] = p.A_grav.y;
		if constexpr (D > 2) AZ_gravs[i] = p.A_grav.z;

		AX_viscs[i] = p.A_visc.x;
		AY_viscs[i] = p.A_visc.y;
		if constexpr (D > 2) AZ_viscs[i] = p.A_visc.z;

		AX_presss[i] = p.A_press.x;
		AY_presss[i] = p.A_press.y;
		if constexpr (D > 2) AZ_presss[i] = p.A_press.z;

		Masss[i] = p.Mass;
		Densitys[i] = p.Density;
		Pressures[i] = p.Pressure;
		Types[i] = p.Type;
		BoundaryPsis[i] = p.BoundaryPsi;
	}
	void SetParticle(size_t i, Particle<D>&& particle) {
		SetParticle(i, particle);
	}
	void SetParticle(size_t i, const ParticleSoA& other, size_t src) {
		Xs[i] = other.Xs[src];
		Ys[i] = other.Ys[src];
		if constexpr (D > 2) Zs[i] = other.Zs[src];

		VXs[i] = other.VXs[src];
		VYs[i] = other.VYs[src];
		if constexpr (D > 2) VZs[i] = other.VZs[src];

		AX_gravs[i] = other.AX_gravs[src];
		AY_gravs[i] = other.AY_gravs[src];
		if constexpr (D > 2) AZ_gravs[i] = other.AZ_gravs[src];

		AX_viscs[i] = other.AX_viscs[src];
		AY_viscs[i] = other.AY_viscs[src];
		if constexpr (D > 2) AZ_viscs[i] = other.AZ_viscs[src];

		AX_presss[i] = other.AX_presss[src];
		AY_presss[i] = other.AY_presss[src];
		if constexpr (D > 2) AZ_presss[i] = other.AZ_presss[src];

		Masss[i] = other.Masss[src];
		Densitys[i] = other.Densitys[src];
		Pressures[i] = other.Pressures[src];
		Types[i] = other.Types[src];
		BoundaryPsis[i] = other.BoundaryPsis[src];
	}


	// -------------------------------------------------
	// Getters
	// -------------------------------------------------

	vec_t Position(size_t i) const {
		vec_t v;
		v.x = Xs[i];
		v.y = Ys[i];
		if constexpr (D > 2) v.z = Zs[i];
		return v;
	}
	float PositionX(size_t i) const { return Xs[i]; }
	float PositionY(size_t i) const { return Ys[i]; }
	float PositionZ(size_t i) const { if constexpr (D > 2) return Zs[i]; else return 0.f; }

	vec_t Velocity(size_t i) const {
		vec_t v;
		v.x = VXs[i];
		v.y = VYs[i];
		if constexpr (D > 2) v.z = VZs[i];
		return v;
	}
	float VelocityX(size_t i) const { return VXs[i]; }
	float VelocityY(size_t i) const { return VYs[i]; }
	float VelocityZ(size_t i) const { if constexpr (D > 2) return VZs[i]; else return 0.f; }

	vec_t AccelerationGravity(size_t i) const {
		vec_t v;
		v.x = AX_gravs[i];
		v.y = AY_gravs[i];
		if constexpr (D > 2) v.z = AZ_gravs[i];
		return v;
	}
	float AccelerationGravityX(size_t i) const { return AX_gravs[i]; }
	float AccelerationGravityY(size_t i) const { return AY_gravs[i]; }
	float AccelerationGravityZ(size_t i) const { if constexpr (D > 2) return AZ_gravs[i]; else return 0.f; }

	vec_t AccelerationViscosity(size_t i) const {
		vec_t v;
		v.x = AX_viscs[i];
		v.y = AY_viscs[i];
		if constexpr (D > 2) v.z = AZ_viscs[i];
		return v;
	}
	float AccelerationViscosityX(size_t i) const { return AX_viscs[i]; }
	float AccelerationViscosityY(size_t i) const { return AY_viscs[i]; }
	float AccelerationViscosityZ(size_t i) const { if constexpr (D > 2) return AZ_viscs[i]; else return 0.f; }

	vec_t AccelerationPressure(size_t i) const {
		vec_t v;
		v.x = AX_presss[i];
		v.y = AY_presss[i];
		if constexpr (D > 2) v.z = AZ_presss[i];
		return v;
	}
	float AccelerationPressureX(size_t i) const { return AX_presss[i]; }
	float AccelerationPressureY(size_t i) const { return AY_presss[i]; }
	float AccelerationPressureZ(size_t i) const { if constexpr (D > 2) return AZ_presss[i]; else return 0.f; }

	float Mass(size_t i) const { return Masss[i]; }
	float Density(size_t i) const { return Densitys[i]; }
	float Pressure(size_t i) const { return Pressures[i]; }
	ParticleType Type(size_t i) const { return Types[i]; }
	float BoundaryPsi(size_t i) const { return BoundaryPsis[i]; }

	// -------------------------------------------------
	// Setters
	// -------------------------------------------------

	void SetPosition(size_t i, const vec_t& p) {
		Xs[i] = p.x;
		Ys[i] = p.y;
		if constexpr (D > 2) Zs[i] = p.z;
	}
	void SetPositionX(size_t i, float v) { Xs[i] = v; }
	void SetPositionY(size_t i, float v) { Ys[i] = v; }
	void SetPositionZ(size_t i, float v) { if constexpr (D > 2) Zs[i] = v; }

	void SetVelocity(size_t i, const vec_t& v) {
		VXs[i] = v.x;
		VYs[i] = v.y;
		if constexpr (D > 2) VZs[i] = v.z;
	}
	void SetVelocityX(size_t i, float v) { VXs[i] = v; }
	void SetVelocityY(size_t i, float v) { VYs[i] = v; }
	void SetVelocityZ(size_t i, float v) { if constexpr (D > 2) VZs[i] = v; }

	void SetAccelerationGravity(size_t i, const vec_t& a) {
		AX_gravs[i] = a.x;
		AY_gravs[i] = a.y;
		if constexpr (D > 2) AZ_gravs[i] = a.z;
	}
	void SetAccelerationGravityX(size_t i, float v) { AX_gravs[i] = v; }
	void SetAccelerationGravityY(size_t i, float v) { AY_gravs[i] = v; }
	void SetAccelerationGravityZ(size_t i, float v) { if constexpr (D > 2) AZ_gravs[i] = v; }

	void SetAccelerationViscosity(size_t i, const vec_t& a) {
		AX_viscs[i] = a.x;
		AY_viscs[i] = a.y;
		if constexpr (D > 2) AZ_viscs[i] = a.z;
	}
	void SetAccelerationViscosityX(size_t i, float v) { AX_viscs[i] = v; }
	void SetAccelerationViscosityY(size_t i, float v) { AY_viscs[i] = v; }
	void SetAccelerationViscosityZ(size_t i, float v) { if constexpr (D > 2) AZ_viscs[i] = v; }

	void SetAccelerationPressure(size_t i, const vec_t& a) {
		AX_presss[i] = a.x;
		AY_presss[i] = a.y;
		if constexpr (D > 2) AZ_presss[i] = a.z;
	}
	void SetAccelerationPressureX(size_t i, float v) { AX_presss[i] = v; }
	void SetAccelerationPressureY(size_t i, float v) { AY_presss[i] = v; }
	void SetAccelerationPressureZ(size_t i, float v) { if constexpr (D > 2) AZ_presss[i] = v; }

	void SetMass(size_t i, float m) { Masss[i] = m; }
	void SetDensity(size_t i, float d) { Densitys[i] = d; }
	void SetPressure(size_t i, float p) { Pressures[i] = p; }
	void SetType(size_t i, ParticleType t) { Types[i] = t; }
	void SetBoundaryPsi(size_t i, float p) { BoundaryPsis[i] = p; }
};


/// Calling Z component functions when D == 2 is undefined behaviour
template <size_t D>
struct ParticleAoS
{
	using vec_t = Particle<D>::vec_t;
	using idx_t = Particle<D>::idx_t;

	std::vector<Particle<D>> ParticlesVector;

	ParticleAoS() = default;
	ParticleAoS(ParticleAoS&& other) : ParticlesVector(std::move(other.ParticlesVector)) {}

	ParticleAoS& operator=(const ParticleAoS& other) { ParticlesVector = other.ParticlesVector; return *this; };
	ParticleAoS& operator=(ParticleAoS&& other) noexcept { ParticlesVector = std::move(other.ParticlesVector); return *this; };


	size_t Size() const { return ParticlesVector.size(); }
	bool Empty() const { return ParticlesVector.empty(); }
	void Resize(size_t size) { ParticlesVector.resize(size); }
	void Clear() { ParticlesVector.clear(); }
	void PushBack(const Particle<D>& particle) { ParticlesVector.push_back(particle); }

	Particle<D> GetParticle(size_t i) const { return ParticlesVector[i]; }
	void GetParticles(std::vector<Particle<D>>& out) const { out = ParticlesVector; }
	void SetParticle(size_t i, const Particle<D>& particle) { ParticlesVector[i] = particle; }
	void SetParticle(size_t i, Particle<D>&& particle) { ParticlesVector[i] = std::move(particle); }
	void SetParticle(size_t i, const ParticleAoS& other, size_t src) { ParticlesVector[i] = other.ParticlesVector[src]; }

	vec_t Position(size_t i)  const { return ParticlesVector[i].Position; }
	float PositionX(size_t i) const { return ParticlesVector[i].Position.x; }
	float PositionY(size_t i) const { return ParticlesVector[i].Position.y; }
	float PositionZ(size_t i) const;

	vec_t Velocity(size_t i)  const { return ParticlesVector[i].Velocity; }
	float VelocityX(size_t i) const { return ParticlesVector[i].Velocity.x; }
	float VelocityY(size_t i) const { return ParticlesVector[i].Velocity.y; }
	float VelocityZ(size_t i) const;

	vec_t AccelerationGravity(size_t i)  const { return ParticlesVector[i].A_grav; }
	float AccelerationGravityX(size_t i) const { return ParticlesVector[i].A_grav.x; }
	float AccelerationGravityY(size_t i) const { return ParticlesVector[i].A_grav.y; }
	float AccelerationGravityZ(size_t i) const;

	vec_t AccelerationViscosity(size_t i)  const { return ParticlesVector[i].A_visc; }
	float AccelerationViscosityX(size_t i) const { return ParticlesVector[i].A_visc.x; }
	float AccelerationViscosityY(size_t i) const { return ParticlesVector[i].A_visc.y; }
	float AccelerationViscosityZ(size_t i) const;

	vec_t AccelerationPressure(size_t i)  const { return ParticlesVector[i].A_press; }
	float AccelerationPressureX(size_t i) const { return ParticlesVector[i].A_press.x; }
	float AccelerationPressureY(size_t i) const { return ParticlesVector[i].A_press.y; }
	float AccelerationPressureZ(size_t i) const;

	float Mass       (size_t i) const { return ParticlesVector[i].Mass; }
	float Density    (size_t i) const { return ParticlesVector[i].Density; }
	float Pressure   (size_t i) const { return ParticlesVector[i].Pressure; }
	ParticleType Type(size_t i) const { return ParticlesVector[i].Type; }
	float BoundaryPsi(size_t i) const { return ParticlesVector[i].BoundaryPsi; }

	void SetPosition(size_t i, const vec_t& pos) { ParticlesVector[i].Position   = pos; }
	void SetPositionX(size_t i, float pos)       { ParticlesVector[i].Position.x = pos; }
	void SetPositionY(size_t i, float pos)       { ParticlesVector[i].Position.y = pos; }
	void SetPositionZ(size_t i, float pos);

	void SetVelocity(size_t i, const vec_t& vel) { ParticlesVector[i].Velocity   = vel; }
	void SetVelocityX(size_t i, float vel)       { ParticlesVector[i].Velocity.x = vel; }
	void SetVelocityY(size_t i, float vel)       { ParticlesVector[i].Velocity.y = vel; }
	void SetVelocityZ(size_t i, float vel);

	void SetAccelerationGravity(size_t i, const vec_t& acc) { ParticlesVector[i].A_grav   = acc; }
	void SetAccelerationGravityX(size_t i, float acc)       { ParticlesVector[i].A_grav.x = acc; }
	void SetAccelerationGravityY(size_t i, float acc)       { ParticlesVector[i].A_grav.y = acc; }
	void SetAccelerationGravityZ(size_t i, float acc);

	void SetAccelerationViscosity(size_t i, const vec_t& acc) { ParticlesVector[i].A_visc   = acc; }
	void SetAccelerationViscosityX(size_t i, float acc)       { ParticlesVector[i].A_visc.x = acc; }
	void SetAccelerationViscosityY(size_t i, float acc)       { ParticlesVector[i].A_visc.y = acc; }
	void SetAccelerationViscosityZ(size_t i, float acc);

	void SetAccelerationPressure(size_t i, const vec_t& acc) { ParticlesVector[i].A_press   = acc; }
	void SetAccelerationPressureX(size_t i, float acc)       { ParticlesVector[i].A_press.x = acc; }
	void SetAccelerationPressureY(size_t i, float acc)       { ParticlesVector[i].A_press.y = acc; }
	void SetAccelerationPressureZ(size_t i, float acc);

	void SetMass(size_t i, float m)           { ParticlesVector[i].Mass        = m;    }
	void SetDensity(size_t i, float d)        { ParticlesVector[i].Density     = d;    }
	void SetPressure(size_t i, float p)       { ParticlesVector[i].Pressure    = p;    }
	void SetType(size_t i, ParticleType type) { ParticlesVector[i].Type        = type; }
	void SetBoundaryPsi(size_t i, float p)    { ParticlesVector[i].BoundaryPsi = p;    }
};


template <>
inline float ParticleAoS<2>::PositionZ(size_t i) const { return 0.0f; }
template <>
inline float ParticleAoS<3>::PositionZ(size_t i) const { return ParticlesVector[i].Position.z; }
template <>
inline float ParticleAoS<2>::VelocityZ(size_t i) const { return 0.0f; }
template <>
inline float ParticleAoS<3>::VelocityZ(size_t i) const { return ParticlesVector[i].Velocity.z; }
template <>
inline float ParticleAoS<2>::AccelerationGravityZ(size_t i) const { return 0.0f; }
template <>
inline float ParticleAoS<3>::AccelerationGravityZ(size_t i) const { return ParticlesVector[i].A_grav.z; }
template <>
inline float ParticleAoS<2>::AccelerationViscosityZ(size_t i) const { return 0.0f; }
template <>
inline float ParticleAoS<3>::AccelerationViscosityZ(size_t i) const { return ParticlesVector[i].A_visc.z; }
template <>
inline float ParticleAoS<2>::AccelerationPressureZ(size_t i) const { return 0.0f; }
template <>
inline float ParticleAoS<3>::AccelerationPressureZ(size_t i) const { return ParticlesVector[i].A_press.z; }

template <>
inline void ParticleAoS<2>::SetPositionZ(size_t i, float pos) {}
template <>
inline void ParticleAoS<3>::SetPositionZ(size_t i, float pos) { ParticlesVector[i].Position.z = pos; }
template <>
inline void ParticleAoS<2>::SetVelocityZ(size_t i, float pos) {}
template <>
inline void ParticleAoS<3>::SetVelocityZ(size_t i, float pos) { ParticlesVector[i].Velocity.z = pos; }
template <>
inline void ParticleAoS<2>::SetAccelerationGravityZ(size_t i, float acc) {}
template <>
inline void ParticleAoS<3>::SetAccelerationGravityZ(size_t i, float acc) { ParticlesVector[i].A_grav.z = acc; }
template <>
inline void ParticleAoS<2>::SetAccelerationViscosityZ(size_t i, float acc) {}
template <>
inline void ParticleAoS<3>::SetAccelerationViscosityZ(size_t i, float acc) { ParticlesVector[i].A_visc.z = acc; }
template <>
inline void ParticleAoS<2>::SetAccelerationPressureZ(size_t i, float acc) {}
template <>
inline void ParticleAoS<3>::SetAccelerationPressureZ(size_t i, float acc) { ParticlesVector[i].A_press.z = acc; }


#ifdef HAS_CUDA
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <device_launch_parameters.h>
#include <cuda.h>

#include "CUDA/Utils.cuh"

/// Calling Z component functions when D == 2 is undefined behaviour
template <size_t D>
struct ParticlesCuda
{
	using vec_t = typename Particle<D>::vec_t;
	using idx_t = typename Particle<D>::idx_t;

	size_t Count = 0;
	float* Xs = nullptr, * Ys = nullptr, * Zs = nullptr;
	float* VXs = nullptr, * VYs = nullptr, * VZs = nullptr;
	float* AX_gravs = nullptr, * AY_gravs = nullptr, * AZ_gravs = nullptr;
	float* AX_viscs = nullptr, * AY_viscs = nullptr, * AZ_viscs = nullptr;
	float* AX_presss = nullptr, * AY_presss = nullptr, * AZ_presss = nullptr;
	float* Masss = nullptr;
	float* Densitys = nullptr;
	float* Pressures = nullptr;
	ParticleType* Types = nullptr;
	float* BoundaryPsis = nullptr;


	ParticlesCuda() = default;

	// -------------------------------------------------
	// Size / Memory
	// -------------------------------------------------

	size_t Size() const { return Count; }
	bool Empty() const { return Count == 0; }
	void Malloc(size_t n)
	{
		// TODO: malloc
		cudaMalloc((void**)&Xs, n * sizeof(float));
		cudaMalloc((void**)&Ys, n * sizeof(float));
		if constexpr (D > 2) cudaMalloc((void**)&Zs, n * sizeof(float));

		cudaMalloc((void**)&VXs, n * sizeof(float));
		cudaMalloc((void**)&VYs, n * sizeof(float));
		if constexpr (D > 2) cudaMalloc((void**)&VZs, n * sizeof(float));

		cudaMalloc((void**)&AX_gravs, n * sizeof(float));
		cudaMalloc((void**)&AY_gravs, n * sizeof(float));
		if constexpr (D > 2) cudaMalloc((void**)&AZ_gravs, n * sizeof(float));

		cudaMalloc((void**)&AX_viscs, n * sizeof(float));
		cudaMalloc((void**)&AY_viscs, n * sizeof(float));
		if constexpr (D > 2) cudaMalloc((void**)&AZ_viscs, n * sizeof(float));

		cudaMalloc((void**)&AX_presss, n * sizeof(float));
		cudaMalloc((void**)&AY_presss, n * sizeof(float));
		if constexpr (D > 2) cudaMalloc((void**)&AZ_presss, n * sizeof(float));

		cudaMalloc((void**)&Masss, n * sizeof(float));
		cudaMalloc((void**)&Densitys, n * sizeof(float));
		cudaMalloc((void**)&Pressures, n * sizeof(float));
		cudaMalloc((void**)&Types, n * sizeof(ParticleType));
		cudaMalloc((void**)&BoundaryPsis, n * sizeof(float));

		Count = n;
	}
	void MemcpyTo(const ParticlesCuda<D>& other) const
	{
		// if we need to resize a non-empty vector: malloc a new one then memcpy to it
	}
	template <ParticleSet<D> Particles>
	void MemcpyTo(Particles& other, u32 stride = 1) const
	{
		if (stride == 1)
		{
			other.Resize(Count);

			cudaMemcpy(other.Xs.data(), Xs, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.Ys.data(), Ys, Count * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2) cudaMemcpy(other.Zs.data(), Zs, Count * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.VXs.data(), VXs, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.VYs.data(), VYs, Count * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2) cudaMemcpy(other.VZs.data(), VZs, Count * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.AX_gravs.data(), AX_gravs, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.AY_gravs.data(), AY_gravs, Count * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2) cudaMemcpy(other.AZ_gravs.data(), AZ_gravs, Count * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.AX_presss.data(), AX_presss, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.AY_presss.data(), AY_presss, Count * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2) cudaMemcpy(other.AZ_presss.data(), AZ_presss, Count * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.AX_viscs.data(), AX_viscs, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.AY_viscs.data(), AY_viscs, Count * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2) cudaMemcpy(other.AZ_viscs.data(), AZ_viscs, Count * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.Masss.data(), Masss, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.Densitys.data(), Densitys, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.Pressures.data(), Pressures, Count * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.Types.data(), Types, Count * sizeof(ParticleType), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.BoundaryPsis.data(), BoundaryPsis, Count * sizeof(float), cudaMemcpyDeviceToHost);
		}
		else
		{
			u64 sample_size = (Count + stride - 1) / stride;

			// Allocate temporary device buffers for sampled data
			float* Xs_sample, * Ys_sample, * Zs_sample = nullptr;
			float* VXs_sample, * VYs_sample, * VZs_sample = nullptr;

			float* AX_gravs_sample, * AY_gravs_sample, * AZ_gravs_sample = nullptr;
			float* AX_presss_sample, * AY_presss_sample, * AZ_presss_sample = nullptr;
			float* AX_viscs_sample, * AY_viscs_sample, * AZ_viscs_sample = nullptr;

			float* Masss_sample, * Densitys_sample, * Pressures_sample;
			float* BoundaryPsis_sample;

			ParticleType* Types_sample;

			// position
			cudaMalloc(&Xs_sample, sample_size * sizeof(float));
			cudaMalloc(&Ys_sample, sample_size * sizeof(float));
			if constexpr (D > 2)
				cudaMalloc(&Zs_sample, sample_size * sizeof(float));

			// velocity
			cudaMalloc(&VXs_sample, sample_size * sizeof(float));
			cudaMalloc(&VYs_sample, sample_size * sizeof(float));
			if constexpr (D > 2)
				cudaMalloc(&VZs_sample, sample_size * sizeof(float));

			// gravity accel
			cudaMalloc(&AX_gravs_sample, sample_size * sizeof(float));
			cudaMalloc(&AY_gravs_sample, sample_size * sizeof(float));
			if constexpr (D > 2)
				cudaMalloc(&AZ_gravs_sample, sample_size * sizeof(float));

			// pressure accel
			cudaMalloc(&AX_presss_sample, sample_size * sizeof(float));
			cudaMalloc(&AY_presss_sample, sample_size * sizeof(float));
			if constexpr (D > 2)
				cudaMalloc(&AZ_presss_sample, sample_size * sizeof(float));

			// viscosity accel
			cudaMalloc(&AX_viscs_sample, sample_size * sizeof(float));
			cudaMalloc(&AY_viscs_sample, sample_size * sizeof(float));
			if constexpr (D > 2)
				cudaMalloc(&AZ_viscs_sample, sample_size * sizeof(float));

			// scalar fields
			cudaMalloc(&Masss_sample, sample_size * sizeof(float));
			cudaMalloc(&Densitys_sample, sample_size * sizeof(float));
			cudaMalloc(&Pressures_sample, sample_size * sizeof(float));
			cudaMalloc(&Types_sample, sample_size * sizeof(ParticleType));
			cudaMalloc(&BoundaryPsis_sample, sample_size * sizeof(float));


			// Launch kernel
			size_t block = 256;
			size_t grid = (sample_size + block - 1) / block;

			SampleStride(grid, block, sample_size, Count, stride,
				Xs, Ys, D > 2 ? Zs : nullptr,
				VXs, VYs, D > 2 ? VZs : nullptr,
				AX_gravs, AY_gravs, D > 2 ? AZ_gravs : nullptr,
				AX_presss, AY_presss, D > 2 ? AZ_presss : nullptr,
				AX_viscs, AY_viscs, D > 2 ? AZ_viscs : nullptr,
				Masss, Densitys, Pressures,
				(u8*)Types, BoundaryPsis,
				Xs_sample, Ys_sample, D > 2 ? Zs_sample : nullptr,
				VXs_sample, VYs_sample, D > 2 ? VZs_sample : nullptr,
				AX_gravs_sample, AY_gravs_sample, D > 2 ? AZ_gravs_sample : nullptr,
				AX_presss_sample, AY_presss_sample, D > 2 ? AZ_presss_sample : nullptr,
				AX_viscs_sample, AY_viscs_sample, D > 2 ? AZ_viscs_sample : nullptr,
				Masss_sample, Densitys_sample, Pressures_sample,
				(u8*)Types_sample, BoundaryPsis_sample,
				D > 2);

			cudaError_t err = cudaDeviceSynchronize();
			if (err != cudaSuccess)
				printf("%s\n", cudaGetErrorString(err));

			cudaError_t cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "SampleKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			}

			// Copy to CPU
			other.Resize(sample_size);

			cudaMemcpy(other.Xs.data(), Xs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.Ys.data(), Ys_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2)
				cudaMemcpy(other.Zs.data(), Zs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.VXs.data(), VXs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.VYs.data(), VYs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2)
				cudaMemcpy(other.VZs.data(), VZs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.AX_gravs.data(), AX_gravs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.AY_gravs.data(), AY_gravs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2)
				cudaMemcpy(other.AZ_gravs.data(), AZ_gravs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.AX_presss.data(), AX_presss_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.AY_presss.data(), AY_presss_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2)
				cudaMemcpy(other.AZ_presss.data(), AZ_presss_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.AX_viscs.data(), AX_viscs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.AY_viscs.data(), AY_viscs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			if constexpr (D > 2)
				cudaMemcpy(other.AZ_viscs.data(), AZ_viscs_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.Masss.data(), Masss_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.Densitys.data(), Densitys_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.Pressures.data(), Pressures_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(other.BoundaryPsis.data(), BoundaryPsis_sample, sample_size * sizeof(float), cudaMemcpyDeviceToHost);

			cudaMemcpy(other.Types.data(), Types_sample, sample_size * sizeof(ParticleType), cudaMemcpyDeviceToHost);
		}
	}
	void MemcpyFrom(const ParticleSoA<D>& other)
	{
		cudaMemcpy(Xs, other.Xs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(Ys, other.Ys.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		if constexpr (D > 2) cudaMemcpy(Zs, other.Zs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(VXs, other.VXs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(VYs, other.VYs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		if constexpr (D > 2) cudaMemcpy(VZs, other.VZs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(AX_gravs, other.AX_gravs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(AY_gravs, other.AY_gravs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		if constexpr (D > 2) cudaMemcpy(AZ_gravs, other.AZ_gravs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(AX_presss, other.AX_presss.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(AY_presss, other.AY_presss.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		if constexpr (D > 2) cudaMemcpy(AZ_presss, other.AZ_presss.data(), Count * sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(AX_viscs, other.AX_viscs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(AY_viscs, other.AY_viscs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		if constexpr (D > 2) cudaMemcpy(AZ_viscs, other.AZ_viscs.data(), Count * sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(Masss, other.Masss.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(Densitys, other.Densitys.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(Pressures, other.Pressures.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(Types, other.Types.data(), Count * sizeof(ParticleType), cudaMemcpyHostToDevice);
		cudaMemcpy(BoundaryPsis, other.BoundaryPsis.data(), Count * sizeof(float), cudaMemcpyHostToDevice);
	}
	void Free()
	{
		if (Empty()) return;

		cudaFree(Xs); cudaFree(Ys); if constexpr (D > 2) cudaFree(Zs);
		cudaFree(VXs); cudaFree(VYs); if constexpr (D > 2) cudaFree(VZs);
		cudaFree(AX_gravs); cudaFree(AY_gravs); if constexpr (D > 2) cudaFree(AZ_gravs);
		cudaFree(AX_viscs); cudaFree(AY_viscs); if constexpr (D > 2) cudaFree(AZ_viscs);
		cudaFree(AX_presss); cudaFree(AY_presss); if constexpr (D > 2) cudaFree(AZ_presss);
		cudaFree(Masss);
		cudaFree(Densitys);
		cudaFree(Pressures);
		cudaFree(Types);
		cudaFree(BoundaryPsis);
	}

	// -------------------------------------------------
	// Conversion
	// -------------------------------------------------
	__device__ Particle<D> DGetParticle(size_t i) const
	{
		Particle<D> p;

		p.Position.x = Xs[i];
		p.Position.y = Ys[i];
		if constexpr (D > 2) p.Position.z = Zs[i];

		p.Velocity.x = VXs[i];
		p.Velocity.y = VYs[i];
		if constexpr (D > 2) p.Velocity.z = VZs[i];

		p.A_grav.x = AX_gravs[i];
		p.A_grav.y = AY_gravs[i];
		if constexpr (D > 2) p.A_grav.z = AZ_gravs[i];

		p.A_visc.x = AX_viscs[i];
		p.A_visc.y = AY_viscs[i];
		if constexpr (D > 2) p.A_visc.z = AZ_viscs[i];

		p.A_press.x = AX_presss[i];
		p.A_press.y = AY_presss[i];
		if constexpr (D > 2) p.A_press.z = AZ_presss[i];

		p.Mass = Masss[i];
		p.Density = Densitys[i];
		p.Pressure = Pressures[i];
		p.Type = Types[i];
		p.BoundaryPsi = BoundaryPsis[i];

		return p;
	}

	__device__ void DSetParticle(size_t i, const Particle<D>& p)
	{
		Xs[i] = p.Position.x;
		Ys[i] = p.Position.y;
		if constexpr (D > 2) Zs[i] = p.Position.z;

		VXs[i] = p.Velocity.x;
		VYs[i] = p.Velocity.y;
		if constexpr (D > 2) VZs[i] = p.Velocity.z;

		AX_gravs[i] = p.A_grav.x;
		AY_gravs[i] = p.A_grav.y;
		if constexpr (D > 2) AZ_gravs[i] = p.A_grav.z;

		AX_viscs[i] = p.A_visc.x;
		AY_viscs[i] = p.A_visc.y;
		if constexpr (D > 2) AZ_viscs[i] = p.A_visc.z;

		AX_presss[i] = p.A_press.x;
		AY_presss[i] = p.A_press.y;
		if constexpr (D > 2) AZ_presss[i] = p.A_press.z;

		Masss[i] = p.Mass;
		Densitys[i] = p.Density;
		Pressures[i] = p.Pressure;
		Types[i] = p.Type;
		BoundaryPsis[i] = p.BoundaryPsi;
	}
	__device__ void DSetParticle(size_t i, Particle<D>&& particle) {
		DSetParticle(i, particle);
	}
	__device__ void DSetParticle(size_t i, const ParticlesCuda& other, size_t src) {
		Xs[i] = other.Xs[src];
		Ys[i] = other.Ys[src];
		if constexpr (D > 2) Zs[i] = other.Zs[src];

		VXs[i] = other.VXs[src];
		VYs[i] = other.VYs[src];
		if constexpr (D > 2) VZs[i] = other.VZs[src];

		AX_gravs[i] = other.AX_gravs[src];
		AY_gravs[i] = other.AY_gravs[src];
		if constexpr (D > 2) AZ_gravs[i] = other.AZ_gravs[src];

		AX_viscs[i] = other.AX_viscs[src];
		AY_viscs[i] = other.AY_viscs[src];
		if constexpr (D > 2) AZ_viscs[i] = other.AZ_viscs[src];

		AX_presss[i] = other.AX_presss[src];
		AY_presss[i] = other.AY_presss[src];
		if constexpr (D > 2) AZ_presss[i] = other.AZ_presss[src];

		Masss[i] = other.Masss[src];
		Densitys[i] = other.Densitys[src];
		Pressures[i] = other.Pressures[src];
		Types[i] = other.Types[src];
		BoundaryPsis[i] = other.BoundaryPsis[src];
	}


	// -------------------------------------------------
	// Getters
	// -------------------------------------------------

	__device__ vec_t DPosition(size_t i) const {
		vec_t v;
		v.x = Xs[i];
		v.y = Ys[i];
		if constexpr (D > 2) v.z = Zs[i];
		return v;
	}
	__device__ float DPositionX(size_t i) const { return Xs[i]; }
	__device__ float DPositionY(size_t i) const { return Ys[i]; }
	__device__ float DPositionZ(size_t i) const { if constexpr (D > 2) return Zs[i]; else return 0.f; }

	__device__ vec_t DVelocity(size_t i) const {
		vec_t v;
		v.x = VXs[i];
		v.y = VYs[i];
		if constexpr (D > 2) v.z = VZs[i];
		return v;
	}
	__device__ float DVelocityX(size_t i) const { return VXs[i]; }
	__device__ float DVelocityY(size_t i) const { return VYs[i]; }
	__device__ float DVelocityZ(size_t i) const { if constexpr (D > 2) return VZs[i]; else return 0.f; }

	__device__ vec_t DAccelerationGravity(size_t i) const {
		vec_t v;
		v.x = AX_gravs[i];
		v.y = AY_gravs[i];
		if constexpr (D > 2) v.z = AZ_gravs[i];
		return v;
	}
	__device__ float DAccelerationGravityX(size_t i) const { return AX_gravs[i]; }
	__device__ float DAccelerationGravityY(size_t i) const { return AY_gravs[i]; }
	__device__ float DAccelerationGravityZ(size_t i) const { if constexpr (D > 2) return AZ_gravs[i]; else return 0.f; }

	__device__ vec_t DAccelerationViscosity(size_t i) const {
		vec_t v;
		v.x = AX_viscs[i];
		v.y = AY_viscs[i];
		if constexpr (D > 2) v.z = AZ_viscs[i];
		return v;
	}
	__device__ float DAccelerationViscosityX(size_t i) const { return AX_viscs[i]; }
	__device__ float DAccelerationViscosityY(size_t i) const { return AY_viscs[i]; }
	__device__ float DAccelerationViscosityZ(size_t i) const { if constexpr (D > 2) return AZ_viscs[i]; else return 0.f; }

	__device__ vec_t DAccelerationPressure(size_t i) const {
		vec_t v;
		v.x = AX_presss[i];
		v.y = AY_presss[i];
		if constexpr (D > 2) v.z = AZ_presss[i];
		return v;
	}
	__device__ float DAccelerationPressureX(size_t i) const { return AX_presss[i]; }
	__device__ float DAccelerationPressureY(size_t i) const { return AY_presss[i]; }
	__device__ float DAccelerationPressureZ(size_t i) const { if constexpr (D > 2) return AZ_presss[i]; else return 0.f; }

	__device__ float DMass(size_t i) const { return Masss[i]; }
	__device__ float DDensity(size_t i) const { return Densitys[i]; }
	__device__ float DPressure(size_t i) const { return Pressures[i]; }
	__device__ ParticleType DType(size_t i) const { return Types[i]; }
	__device__ float DBoundaryPsi(size_t i) const { return BoundaryPsis[i]; }

	// -------------------------------------------------
	// Setters
	// -------------------------------------------------

	__device__ void DSetPosition(size_t i, const vec_t& p) {
		Xs[i] = p.x;
		Ys[i] = p.y;
		if constexpr (D > 2) Zs[i] = p.z;
	}
	__device__ void DSetPositionX(size_t i, float v) { Xs[i] = v; }
	__device__ void DSetPositionY(size_t i, float v) { Ys[i] = v; }
	__device__ void DSetPositionZ(size_t i, float v) { if constexpr (D > 2) Zs[i] = v; }

	__device__ void DSetVelocity(size_t i, const vec_t& v) {
		VXs[i] = v.x;
		VYs[i] = v.y;
		if constexpr (D > 2) VZs[i] = v.z;
	}
	__device__ void DSetVelocityX(size_t i, float v) { VXs[i] = v; }
	__device__ void DSetVelocityY(size_t i, float v) { VYs[i] = v; }
	__device__ void DSetVelocityZ(size_t i, float v) { if constexpr (D > 2) VZs[i] = v; }

	__device__ void DSetAccelerationGravity(size_t i, const vec_t& a) {
		AX_gravs[i] = a.x;
		AY_gravs[i] = a.y;
		if constexpr (D > 2) AZ_gravs[i] = a.z;
	}
	__device__ void DSetAccelerationGravityX(size_t i, float v) { AX_gravs[i] = v; }
	__device__ void DSetAccelerationGravityY(size_t i, float v) { AY_gravs[i] = v; }
	__device__ void DSetAccelerationGravityZ(size_t i, float v) { if constexpr (D > 2) AZ_gravs[i] = v; }

	__device__ void DSetAccelerationViscosity(size_t i, const vec_t& a) {
		AX_viscs[i] = a.x;
		AY_viscs[i] = a.y;
		if constexpr (D > 2) AZ_viscs[i] = a.z;
	}
	__device__ void DSetAccelerationViscosityX(size_t i, float v) { AX_viscs[i] = v; }
	__device__ void DSetAccelerationViscosityY(size_t i, float v) { AY_viscs[i] = v; }
	__device__ void DSetAccelerationViscosityZ(size_t i, float v) { if constexpr (D > 2) AZ_viscs[i] = v; }

	__device__ void DSetAccelerationPressure(size_t i, const vec_t& a) {
		AX_presss[i] = a.x;
		AY_presss[i] = a.y;
		if constexpr (D > 2) AZ_presss[i] = a.z;
	}
	__device__ void DSetAccelerationPressureX(size_t i, float v) { AX_presss[i] = v; }
	__device__ void DSetAccelerationPressureY(size_t i, float v) { AY_presss[i] = v; }
	__device__ void DSetAccelerationPressureZ(size_t i, float v) { if constexpr (D > 2) AZ_presss[i] = v; }

	__device__ void DSetMass(size_t i, float m) { Masss[i] = m; }
	__device__ void DSetDensity(size_t i, float d) { Densitys[i] = d; }
	__device__ void DSetPressure(size_t i, float p) { Pressures[i] = p; }
	__device__ void DSetType(size_t i, ParticleType t) { Types[i] = t; }
	__device__ void DSetBoundaryPsi(size_t i, float p) { BoundaryPsis[i] = p; }
};
#endif