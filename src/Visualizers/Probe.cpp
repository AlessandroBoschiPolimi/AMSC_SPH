#include "Probe.hpp"

int Probe::NEXT_ID = 0;

std::optional<Probe> Probe::ParseOne(const std::string& str)
{
	std::stringstream ss(str);
	std::string token;
	Probe p;

	// Extract the name
	if (!std::getline(ss, token, ':'))
		return std::nullopt;
	p.Name = token;

	auto parseFloat = [](const std::string& str, float& out) -> bool {
			try {
				out = std::stof(str);
				return true;
			}
			catch (const std::invalid_argument&) {
				return false;
			}
			catch (const std::out_of_range&) {
				return false;
			}
		};

	// Extract the coords
	if (!std::getline(ss, token, ':') || !parseFloat(token, p.TL.x)) { return std::nullopt; }
	if (!std::getline(ss, token, ':') || !parseFloat(token, p.TL.y)) { return std::nullopt; }
	if (!std::getline(ss, token, ':') || !parseFloat(token, p.BR.x)) { return std::nullopt; }
	if (!std::getline(ss, token, ';') || !parseFloat(token, p.BR.y)) { return std::nullopt; }

	return p;
}
std::optional<std::vector<Probe>> Probe::ParseMultiple(const std::string& str)
{
	std::vector<Probe> ps;
	if (str.empty())
		return ps;

	std::stringstream ss(str);
	std::string token;

	// Get one probe data
	if (!std::getline(ss, token, ';'))
		return std::nullopt;

	do
	{
		std::optional<Probe> p = Probe::ParseOne(token);
		if (p.has_value())
			ps.push_back(p.value());
		else
			return std::nullopt;
	} while (std::getline(ss, token, ';'));

	return ps;
}

std::string Probe::ToString() const
{
	std::stringstream ss;

	ss << Name << ':';
	ss << TL.x << ':' << TL.y << ':';
	ss << BR.x << ':' << BR.y;

	return ss.str();
}