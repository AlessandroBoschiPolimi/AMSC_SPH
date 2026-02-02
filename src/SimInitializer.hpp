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


template <size_t D>
class BoxInitializer;

template <>
class BoxInitializer<2> : public SimInitializer<2>
{
public:
	void Init(std::vector<Particle<2>>& out, float h) const override;
};

inline void BoxInitializer<2>::Init(std::vector<Particle<2>>& out, float h) const
{
	using vec_t = Particle<2>::vec_t;

	// Setup for "Dam breaking"
	int maxx = 100, maxy = 100;
	vec_t vertical_direction{ 0, -1 };
	vec_t zero_direction{ 0, 0 };
	for (float x = 2; x < maxx / 2; x += 0.5)
	{
		for (float y = 28; y < 78; y += 0.5)
		{
			Particle<2> p;
			p.Type = FLUID;
			p.Position.y = float(y) / maxy;
			p.Position.x = float(x) / maxx;
			p.Velocity = zero_direction;
			// Arbitrary Mass to match with other constants
			p.Mass = 0.5f;
			p.A_grav = G_CONSTANT * vertical_direction;
			out.push_back(p);
		}
	}

	// Create multiple walls of the box to avoid leaks
	// TODO: Add walls options to properties
	for (int i = 0; i < 5; i++)
	{
		BuildWallAlongX(out, 0.25 - h / 400 * i, 0   , 1.0, 0.0025);
		BuildWallAlongX(out, 0.8  + h / 400 * i, 0   , 1.0, 0.0025);
		BuildWallAlongY(out, 0    + h / 400 * i, 0.25, 0.8, 0.0025);
		BuildWallAlongY(out, 1.0  - h / 400 * i, 0.25, 0.8, 0.0025);
	}
}



template <size_t D>
class TrayInitializer;

template <>
class TrayInitializer<3> : public SimInitializer<3>
{
public:
	void Init(std::vector<Particle<3>>& out, float h) const override;
};

inline void TrayInitializer<3>::Init(std::vector<Particle<3>>& out, float h) const
{
	using vec_t = Particle<3>::vec_t;

	// Setup for "Dam breaking"
	{
		int maxx = 100, maxy = 100, maxz = 100;
		vec_t vertical_direction{ 0, 0, -1 };
		vec_t zero_direction{ 0, 0, 0 };
		for (float x = 5; x < maxx / 5; x += 0.5)
		{
			for (float y = 5; y < maxy / 5; y += 0.5)
			{
				for (float z = 5; z < maxz / 5; z += 0.5)
				{
					Particle<3> p;
					p.Type = FLUID;
					p.Position.x = float(x) / maxx;
					p.Position.y = float(y) / maxy;
					p.Position.z = float(z) / maxz;
					p.Velocity = zero_direction;
					// Arbitrary Mass to match with other constants
					p.Mass = 0.5f;
					p.A_grav = G_CONSTANT * vertical_direction;
					out.push_back(p);
				}
			}
		}
	}

	// Create multiple walls of the box to avoid leaks
	// TODO: Add walls options to properties
	{
		float maxx = 1.0 / 4, maxy = 1.0 / 4, maxz = 0.1;
		for (int i = 0; i < 5; i++)
		{
			BuildWallAlongXY(out, 0 + h / 400 * i, 0, maxx, 0, maxy, 0.0025);
			BuildWallAlongYZ(out, 0 + h / 400 * i, 0, maxx, 0, maxz, 0.0025);
			BuildWallAlongYZ(out, maxx - h / 400 * i, 0, maxx, 0, maxz, 0.0025);
			BuildWallAlongXZ(out, 0 + h / 400 * i, 0, maxy, 0, maxz, 0.0025);
			BuildWallAlongXZ(out, maxy - h / 400 * i, 0, maxy, 0, maxz, 0.0025);
		}
	}
}





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