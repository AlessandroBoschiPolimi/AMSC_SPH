#pragma once
#include "Utility.hpp"

enum ParticleType : u8
{
	FLUID,
	SOLID
};

/*
AoS   : [{x, y, z, vx, vy, vz, ...}]             Array of particles
Hybrid: {[{x, y, z}], [{vx, vy, vz}], ...}       One array per coord
SoA   : {[x], [y], [z], [vx], [vy], [vz], ...}   One array per component
*/


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

template <typename P, size_t D>
concept ParticleSet = std::movable<P> && std::default_initializable<P> &&
requires(const P cp, P p, P other, size_t i, const Particle<D>& particle, Particle<D>&& particle_move)
{
	{ cp.Size()				} -> std::convertible_to<size_t>;
	{ cp.Empty()			} -> std::convertible_to<bool>;
	{ p.Resize(size_t{})    };
	{ p.Clear()				};
	{ p.PushBack(particle)	};

	{ cp.GetParticle(i)					} -> std::same_as<Particle<D>>;
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

//template <size_t D>
//struct Particles
//{
//	struct Particle
//	{
//		vec_t Position;
//		vec_t Velocity;
//		vec_t A_grav;
//		vec_t A_visc;
//		vec_t A_press;
//		float Mass;
//		float Density;
//		float Pressure;
//		ParticleType Type;
//		float BoundaryPsi;
//	};
//
//	virtual bool Empty() = 0;
//	virtual size_t Size() = 0;
//
//	virtual float PositionX(idx_t idx) = 0;
//	virtual float PositionY(idx_t idx) = 0;
//	virtual float PositionZ(idx_t idx) = 0; // TODO: only for 3D, or make it return garbage
//	virtual vec_t Position (idx_t idx) = 0;
//	virtual ParticleType Type(idx_t idx) = 0;
//};


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
};


template <size_t D>
struct ParticleSoA
{
	using vec_t = Particle<D>::vec_t;
	using idx_t = Particle<D>::idx_t;

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
};


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

	float Mass(size_t i) const { return ParticlesVector[i].Mass; }
	float Density(size_t i) const { return ParticlesVector[i].Density; }
	float Pressure(size_t i) const { return ParticlesVector[i].Pressure; }
	ParticleType Type(size_t i) const { return ParticlesVector[i].Type; }
	float BoundaryPsi(size_t i) const { return ParticlesVector[i].BoundaryPsi; }

	void SetPosition(size_t i, const vec_t& pos) { ParticlesVector[i].Position = pos; }
	void SetPositionX(size_t i, float pos) { ParticlesVector[i].Position.x = pos; }
	void SetPositionY(size_t i, float pos) { ParticlesVector[i].Position.y = pos; }
	void SetPositionZ(size_t i, float pos);

	void SetVelocity(size_t i, const vec_t& vel) { ParticlesVector[i].Velocity = vel; }
	void SetVelocityX(size_t i, float vel) { ParticlesVector[i].Velocity.x = vel; }
	void SetVelocityY(size_t i, float vel) { ParticlesVector[i].Velocity.y = vel; }
	void SetVelocityZ(size_t i, float vel);

	void SetAccelerationGravity(size_t i, const vec_t& acc) { ParticlesVector[i].A_grav = acc; }
	void SetAccelerationGravityX(size_t i, float acc) { ParticlesVector[i].A_grav.x = acc; }
	void SetAccelerationGravityY(size_t i, float acc) { ParticlesVector[i].A_grav.y = acc; }
	void SetAccelerationGravityZ(size_t i, float acc);

	void SetAccelerationViscosity(size_t i, const vec_t& acc) { ParticlesVector[i].A_visc = acc; }
	void SetAccelerationViscosityX(size_t i, float acc) { ParticlesVector[i].A_visc.x = acc; }
	void SetAccelerationViscosityY(size_t i, float acc) { ParticlesVector[i].A_visc.y = acc; }
	void SetAccelerationViscosityZ(size_t i, float acc);

	void SetAccelerationPressure(size_t i, const vec_t& acc) { ParticlesVector[i].A_press = acc; }
	void SetAccelerationPressureX(size_t i, float acc) { ParticlesVector[i].A_press.x = acc; }
	void SetAccelerationPressureY(size_t i, float acc) { ParticlesVector[i].A_press.y = acc; }
	void SetAccelerationPressureZ(size_t i, float acc);

	void SetMass(size_t i, float m) { ParticlesVector[i].Mass = m; }
	void SetDensity(size_t i, float d) { ParticlesVector[i].Density = d; }
	void SetPressure(size_t i, float p) { ParticlesVector[i].Pressure = p; }
	void SetType(size_t i, ParticleType type) { ParticlesVector[i].Type = type; }
	void SetBoundaryPsi(size_t i, float p) { ParticlesVector[i].BoundaryPsi = p; }
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