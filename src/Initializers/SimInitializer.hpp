#pragma once
#include <memory>

#include "Objects/Object.hpp"
#include "Objects/Sink.hpp"
#include "Objects/Source.hpp"
#include "Particle.hpp"

template <size_t D, ParticleSet<D> Particles>
class SimInitializer;

template <ParticleSet<2> Particles>
class SimInitializer<2, Particles>
{
public:
	using Objects = std::vector<uptr<Object<2, Particles>>>;

	SimInitializer() = default;
	virtual ~SimInitializer() = default;

	virtual void Init(Particles& out, float h, Objects& obj) const = 0;

	void BuildWallAlongX(Particles& out, float y, float minx, float maxx, float delta) const;
	void BuildWallAlongY(Particles& out, float x, float miny, float maxy, float delta) const;
	void BuildBox(Particles& out, coord<float, 2> left_bottom_corner, float xsize, float ysize, float delta) const;
	void BuildCircle(Particles& out, coord<float, 2> centre, float radius, float delta) const;
	void AddSink(Objects& obj, const coord<float, 2> A, const coord<float, 2> B, Particles& out, const bool is_right) const;
	void AddSource(Objects& obj, const coord<float, 2> A, const coord<float, 2> B, Particles& out, const bool is_right, const float v, const int parts, const int fp_count, const float mass, const bool is_gravity, const int max_frame = 0) const;

	/// Returns { min-xy, max-xy }
	virtual std::pair<coord<float, 2>, coord<float, 2>> GetDomain() const = 0;
};

template <ParticleSet<3> Particles>
class SimInitializer<3, Particles>
{
public:
	using Objects = std::vector<uptr<Object<3, Particles>>>;

	SimInitializer() = default;
	virtual ~SimInitializer() = default;

	virtual void Init(Particles& out, float h, Objects& obj) const = 0;

	void BuildWallAlongXZ(Particles& out, float y, float minx, float maxx, float minz, float maxz, float delta) const;
	void BuildWallAlongYZ(Particles& out, float x, float miny, float maxy, float minz, float maxz, float delta) const;
	void BuildWallAlongXY(Particles& out, float z, float minx, float maxx, float miny, float maxy, float delta) const;

	/// Returns { min-xyz, max-xyz }
	virtual std::pair<coord<float, 3>, coord<float, 3>> GetDomain() const = 0;
};







template <ParticleSet<2> Particles>
inline void SimInitializer<2, Particles>::BuildWallAlongX(Particles& out, float y, float minx, float maxx, float delta) const
{
	for (float x = minx; x <= maxx; x += delta)
	{
		Particle<2> p;
		p.Type = SOLID;
		p.Position.y = y;
		p.Position.x = x;
		p.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.PushBack(p);
	}
}
template <ParticleSet<2> Particles>
inline void SimInitializer<2, Particles>::BuildWallAlongY(Particles& out, float x, float miny, float maxy, float delta) const
{
	for (float y = miny; y <= maxy; y += delta)
	{
		Particle<2> p;
		p.Type = SOLID;
		p.Position.y = y;
		p.Position.x = x;
		p.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.PushBack(p);
	}
}
template <ParticleSet<2> Particles>
inline void SimInitializer<2, Particles>::BuildBox(Particles& out, coord<float, 2> left_bottom_corner, float xsize, float ysize, float delta) const
{
	for (float x = 0; x <= xsize; x+= delta)
	{
		Particle<2> p_bot;
		p_bot.Type = SOLID;
		p_bot.Position.y = left_bottom_corner.y;
		p_bot.Position.x = left_bottom_corner.x + x;
		p_bot.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.PushBack(p_bot);

		Particle<2> p_top;
		p_top.Type = SOLID;
		p_top.Position.y = left_bottom_corner.y + ysize;
		p_top.Position.x = left_bottom_corner.x + x;
		p_top.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.PushBack(p_top);
	}
	for (float y = delta; y < ysize; y+= delta)
	{
		Particle<2> p_bot;
		p_bot.Type = SOLID;
		p_bot.Position.y = left_bottom_corner.y + y;
		p_bot.Position.x = left_bottom_corner.x;
		p_bot.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.PushBack(p_bot);

		Particle<2> p_top;
		p_top.Type = SOLID;
		p_top.Position.y = left_bottom_corner.y + y;
		p_top.Position.x = left_bottom_corner.x + xsize;
		p_top.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.PushBack(p_top);
	}
}
template <ParticleSet<2> Particles>
inline void SimInitializer<2, Particles>::BuildCircle(Particles& out, coord<float, 2> centre, float radius, float delta) const
{
	for (float psi = 0.0f; psi <= 2.0f * std::numbers::pi; psi += delta)
	{
		Particle<2> p;
		p.Type = SOLID;
		p.Position.y = centre.y + std::sin(psi) * radius;
		p.Position.x = centre.x + std::cos(psi) * radius;
		p.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.PushBack(p);
	}
}

template <ParticleSet<2> Particles>
inline void SimInitializer<2, Particles>::AddSink(Objects& obj, const coord<float, 2> A, const coord<float, 2> B, Particles& out, const bool is_right) const
{
	obj.push_back(std::make_unique<Sink<2, Particles>>(A, B, out, is_right));
}
template <ParticleSet<2> Particles>
inline void SimInitializer<2, Particles>::AddSource(Objects& obj, const coord<float, 2> A, const coord<float, 2> B, Particles& out, const bool is_right, const float v, const int parts, const int fp_count, const float mass, const bool is_gravity, const int max_frame) const
{
	auto o = std::make_unique<Source<2, Particles>>(A, B, out, is_right);
	o->SetParams(v, parts, fp_count, mass, is_gravity, max_frame);
	obj.push_back(std::move(o));
}




template <ParticleSet<3> Particles>
inline void SimInitializer<3, Particles>::BuildWallAlongXY(Particles& out, float z, float minx, float maxx, float miny, float maxy, float delta) const
{
	for (float x = minx; x <= maxx; x += delta)
	{
		for (float y = miny; y <= maxy; y += delta)
		{
			Particle<3> p;
			p.Type = SOLID;
			p.Position.z = z;
			p.Position.y = y;
			p.Position.x = x;
			p.Velocity = Particle<3>::vec_t{ 0.0f, 0.0f, 0.0f };
			out.PushBack(p);
		}
	}
}
template <ParticleSet<3> Particles>
inline void SimInitializer<3, Particles>::BuildWallAlongYZ(Particles& out, float x, float miny, float maxy, float minz, float maxz, float delta) const
{
	for (float z = minz; z <= maxz; z += delta)
	{
		for (float y = miny; y <= maxy; y += delta)
		{
			Particle<3> p;
			p.Type = SOLID;
			p.Position.z = z;
			p.Position.y = y;
			p.Position.x = x;
			p.Velocity = Particle<3>::vec_t{ 0.0f, 0.0f, 0.0f };
			out.PushBack(p);
		}
	}
}
template <ParticleSet<3> Particles>
inline void SimInitializer<3, Particles>::BuildWallAlongXZ(Particles& out, float y, float minx, float maxx, float minz, float maxz, float delta) const
{
	for (float x = minx; x <= maxx; x += delta)
	{
		for (float z = minz; z <= maxz; z += delta)
		{
			Particle<3> p;
			p.Type = SOLID;
			p.Position.z = z;
			p.Position.y = y;
			p.Position.x = x;
			p.Velocity = Particle<3>::vec_t{ 0.0f, 0.0f, 0.0f };
			out.PushBack(p);
		}
	}
}
