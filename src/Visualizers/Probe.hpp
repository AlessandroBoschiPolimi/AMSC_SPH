#pragma once
#include "Utility.hpp"

struct Probe
{
	static int NEXT_ID;
	static const int INVALID_ID = -1;

	const int ID = INVALID_ID;
	std::string Name;
	bool Selected = false;
	coord<float, 2> TL, BR;

	Probe() : ID(NEXT_ID++) {}

	Probe(const Probe& other)
		: ID(NEXT_ID++)
		, Name(other.Name)
		, TL(other.TL)
		, BR(other.BR)
		, Selected(other.Selected)
	{}
	Probe(Probe&& other)
		: ID(other.ID)
		, Name(other.Name)
		, TL(other.TL)
		, BR(other.BR)
		, Selected(other.Selected)
	{}

	void operator=(const Probe& other)
	{
		Name = other.Name;
		TL = other.TL;
		BR = other.BR;
		Selected = other.Selected;
	}

	static std::optional<Probe> ParseOne(const std::string& str);
	static std::optional<std::vector<Probe>> ParseMultiple(const std::string& str);

	std::string ToString() const;
};