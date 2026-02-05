#pragma once
#include "Particle.hpp"

template <size_t D>
class SimInitializer;

template <>
class SimInitializer<2>
{
public:
	SimInitializer() = default;
	virtual ~SimInitializer() = default;

	virtual void Init(std::vector<Particle<2>>& out, float h) const = 0;

	void BuildWallAlongX(std::vector<Particle<2>>& out, float y, float minx, float maxx, float delta) const;
	void BuildWallAlongY(std::vector<Particle<2>>& out, float x, float miny, float maxy, float delta) const;
};

template <>
class SimInitializer<3>
{
public:
	SimInitializer() = default;
	virtual ~SimInitializer() = default;

	virtual void Init(std::vector<Particle<3>>& out, float h) const = 0;

	void BuildWallAlongXZ(std::vector<Particle<3>>& out, float y, float minx, float maxx, float minz, float maxz, float delta) const;
	void BuildWallAlongYZ(std::vector<Particle<3>>& out, float x, float miny, float maxy, float minz, float maxz, float delta) const;
	void BuildWallAlongXY(std::vector<Particle<3>>& out, float z, float minx, float maxx, float miny, float maxy, float delta) const;
};







inline void SimInitializer<2>::BuildWallAlongX(std::vector<Particle<2>>& out, float y, float minx, float maxx, float delta) const
{
	for (float x = minx; x <= maxx; x += delta)
	{
		Particle<2> p;
		p.Type = SOLID;
		p.Position.y = y;
		p.Position.x = x;
		p.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.push_back(p);
	}
}
inline void SimInitializer<2>::BuildWallAlongY(std::vector<Particle<2>>& out, float x, float miny, float maxy, float delta) const
{
	for (float y = miny; y <= maxy; y += delta)
	{
		Particle<2> p;
		p.Type = SOLID;
		p.Position.y = y;
		p.Position.x = x;
		p.Velocity = Particle<2>::vec_t{ 0.0f, 0.0f };
		out.push_back(p);
	}
}


inline void SimInitializer<3>::BuildWallAlongXY(std::vector<Particle<3>>& out, float z, float minx, float maxx, float miny, float maxy, float delta) const
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
			out.push_back(p);
		}
	}
}
inline void SimInitializer<3>::BuildWallAlongYZ(std::vector<Particle<3>>& out, float x, float miny, float maxy, float minz, float maxz, float delta) const
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
			out.push_back(p);
		}
	}
}
inline void SimInitializer<3>::BuildWallAlongXZ(std::vector<Particle<3>>& out, float y, float minx, float maxx, float minz, float maxz, float delta) const
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
			out.push_back(p);
		}
	}
}