#include "Utility.hpp"


void print_time(const stdc::nanoseconds& time, std::ostream& out)
{
	double value = to<double>(time.count());
	const char* unit = "ns";

	if (time >= 1s) {
		value /= 1'000'000'000.0;
		unit = "s";
	}
	else if (time >= 1ms) {
		value /= 1'000'000.0;
		unit = "ms";
	}
	else if (time >= 1us) {
		value /= 1'000.0;
		unit = "us";
	}

#ifdef HAS_CPP20
	out << std::format("{:.2f}{}", value, unit);
#else
	std::ios oldState(nullptr);
	oldState.copyfmt(out);

	out << std::fixed << std::setprecision(2) << value << unit;

	out.copyfmt(oldState);
#endif
}


#ifdef HAS_CPP17
std::string ReadFile(const fs::path& file)
{
	std::ifstream in(file);
	std::string content = std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
	in.close();
	return content;
}
#endif
