/*
Per me si va nella citta' bacata,
per me si va nell'etterno ciclo,
per me si va tra i comportamenti indefiniti.

Giustizia mosse il mio alto fattore; [Bjarne Stroustrup]
fecemi la divina podestate,
la somma sapienza e 'l primo amore;

Dinanzi a me non fuor cose create
se non etterne, e io etterno duro.
Lasciate ogne ';', voi ch'intrate.

-- Iscrizione presente nell'entry point di ogni compiler per cpp
-- Dante Programmieri, Divina Esecuzione, Compilatore, Canto III
*/

#pragma once

#if defined(_MSVC_LANG) && _MSVC_LANG > __cplusplus
#define CPP_VERSION _MSVC_LANG
#else
#define CPP_VERSION __cplusplus
#endif

#if CPP_VERSION >= 202302L
	#define HAS_CPP23
#endif
#if CPP_VERSION >= 202002L
	#define HAS_CPP20
#endif
#if CPP_VERSION >= 201703L
	#define HAS_CPP17
#endif
#if CPP_VERSION >= 201402L
	#define HAS_CPP14
#else
	#error "Bro update your computer"
#endif

#ifdef HAS_CPP20
	#define _requires(...) requires __VA_ARGS__
#else
	#define _requires(...)
#endif


#include <inttypes.h>
#include <numbers>

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <stack>
#include <queue>

#include <functional>
#include <cmath>

#include <memory>
#include <condition_variable>
#include <semaphore>

#include <regex>
#include <chrono>

#ifdef HAS_CPP20
#include <ranges>
#endif
#ifdef HAS_CPP17
#include <filesystem>
#include <optional>
#include <variant>
#endif

#ifdef _DEBUG
#define LOG(x)  std::cout << x
#define LOGE(x) std::cerr << x
#else
#define LOG(x)
#define LOGE(x)
#endif

// ################################################################## ALIASES ##################################################################
#ifdef HAS_CPP20
namespace stdv = std::views;
namespace stdr = std::ranges;
#endif
#ifdef HAS_CPP17
namespace fs   = std::filesystem;
#endif
namespace stdc = std::chrono;
using namespace std::string_literals;
using namespace std::chrono_literals;

using u8  = uint8_t ;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i8  =  int8_t ;
using i16 =  int16_t;
using i32 =  int32_t;
using i64 =  int64_t;
using f32 =  float  ;
using f64 =  double ;

template <typename K, typename V, class H = std::hash<K>>
using hmap = std::unordered_map<K, V, H>;
template <typename V>
using hset = std::unordered_set<V>;

using stdclock = stdc::high_resolution_clock;

template <typename T>
using uptr = std::unique_ptr<T>;
template <typename T>
using sptr = std::shared_ptr<T>;

using callable = std::function<void(void)>;
template <typename T>
using consumer = std::function<void(T)>;
template <typename T>
using producer = std::function<T(void)>;
// ################################################################## ALIASES ##################################################################

static constexpr float G_CONSTANT = 9.81f;

// CASTING
template <typename T, typename X>
inline constexpr T to(const X& x) { return static_cast<T>(x); }
template <typename T, typename X, typename _X>
inline constexpr T to(const stdc::duration<_X, X>& x) { return stdc::duration_cast<T>(x); }


// ################################################################## COORD ##################################################################
template <typename T, size_t D>
struct coord;

template <typename T>
struct coord<T, 2> {
	using type = T;
	T x, y;
};
template <typename T>
struct coord<T, 3> {
	using type = T;
	T x, y, z;
};


template <typename T>
inline coord<T, 2> operator/(const coord<T, 2>& a, const T C) {
	return { a.x / C,  a.y / C };
}
template <typename T>
inline coord<T, 3> operator/(const coord<T, 3>& a, const T C) {
	return { a.x / C,  a.y / C,  a.z / C };
}
template <typename T>
inline coord<T, 2> operator*(const T C, const coord<T, 2>& a) {
	return { C * a.x, C * a.y };
}
template <typename T>
inline coord<T, 3> operator*(const T C, const coord<T, 3>& a) {
	return { C * a.x, C * a.y, C * a.z };
}

template <typename T>
inline coord<T, 2> operator+(const coord<T, 2>& a, const coord<T, 2>& b) {
	return { a.x + b.x, a.y + b.y };
}
template <typename T>
inline coord<T, 3> operator+(const coord<T, 3>& a, const coord<T, 3>& b) {
	return { a.x + b.x, a.y + b.y, a.z + b.z };
}
template <typename T>
inline coord<T, 2> operator-(const coord<T, 2>& a, const coord<T, 2>& b) {
	return { a.x - b.x, a.y - b.y };
}
template <typename T>
inline coord<T, 3> operator-(const coord<T, 3>& a, const coord<T, 3>& b) {
	return { a.x - b.x, a.y - b.y, a.z - b.z };
}

template <typename T, size_t D>
inline coord<T, D>& operator+=(coord<T, D>& a, const coord<T, D>& b) {
	a = a + b;
	return a;
}
template <typename T, size_t D>
inline coord<T, D>& operator-=(coord<T, D>& a, const coord<T, D>& b) {
	a = a - b;
	return a;
}

template <typename T>
inline bool operator==(const coord<T, 2>& a, const coord<T, 2>& b) {
	return a.x == b.x && a.y == b.y;
}
template <typename T>
inline bool operator==(const coord<T, 3>& a, const coord<T, 3>& b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

template <typename T, typename R = T>
inline R Dot(const coord<T, 2>& a, const coord<T, 2>& b) {
	return to<R>(a.x) * to<R>(b.x) + to<R>(a.y) * to<R>(b.y);
}
template <typename T, typename R = T>
inline R Dot(const coord<T, 3>& a, const coord<T, 3>& b) {
	return to<R>(a.x) * to<R>(b.x) + to<R>(a.y) * to<R>(b.y) + to<R>(a.z) * to<R>(b.z);
}
template <typename T, size_t D, typename R = T>
inline R Norm(const coord<T, D>& a) {
	return std::sqrt(Dot(a, a));
}


template <typename T, typename R>
inline constexpr T to(const coord<R, 2>& c) {
	return { to<typename T::type>(c.x), to<typename T::type>(c.y) };
}
template <typename T, typename R>
inline constexpr T to(const coord<R, 3>& c) {
	return { to<typename T::type>(c.x), to<typename T::type>(c.y), to<typename T::type>(c.z) };
}


template <size_t D>
struct CoordIntHash;

template <>
struct CoordIntHash<2> {
	size_t operator()(const coord<int, 2>& p) const {
		size_t h1 = std::hash<int>{}(p.x);
		size_t h2 = std::hash<int>{}(p.y);

		return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
	}
};
template <>
struct CoordIntHash<3> {
	size_t operator()(const coord<int, 3>& p) const {
		return ((to<size_t>(p.x) * 73856093) ^
				(to<size_t>(p.y) * 19349663) ^
				(to<size_t>(p.z) * 83492791));

		size_t h1 = std::hash<int>{}(p.x);
		size_t h2 = std::hash<int>{}(p.y);
		size_t h3 = std::hash<int>{}(p.z);

		size_t seed = h1;
		seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
		seed ^= h3 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
		return seed;
	}
};






template <typename T, size_t D>
struct ParticlesVec;

template <typename T>
struct ParticlesVec<T, 2> {
	std::vector<T> x, y;

	inline coord<T, 2>& operator[](const size_t i) {
		return { x[i], y[i] };
	}
};
template <typename T>
struct ParticlesVec<T, 3> {
	std::vector<T> x, y, z;

	inline coord<T, 3>& operator[](const size_t i) {
		return { x[i], y[i], z[i] };
	}
};
// ################################################################## COORD ##################################################################




// CONCEPTS
#ifdef HAS_CPP17
template <typename T, typename NonVariant>
struct variant_contains : std::false_type {};
template <typename T, typename... Types>
struct variant_contains<T, std::variant<Types...>> : std::disjunction<std::is_same<T, Types>...> {};
#endif
#ifdef HAS_CPP20
template <typename T, typename Variant>
concept in_variant = variant_contains<T, Variant>::value;
#endif

#ifdef HAS_CPP17
template <typename T, typename... Args>
constexpr bool is_one_of = (std::is_same_v<T, Args> || ...);
template <auto Value, auto... Accepted>
constexpr bool val_is_one_of = ((Value == Accepted) || ...);
#endif

template <typename>
struct is_std_function : std::false_type {};
template <typename R, typename... Args>
struct is_std_function<std::function<R(Args...)>> : std::true_type {};
#ifdef HAS_CPP20
template <typename F>
concept StdFunction = is_std_function<F>::value;
#endif

template <typename>
struct c_function_pointer;
template <typename R, typename... Args>
struct c_function_pointer<std::function<R(Args...)>> {
	using type = R(*)(Args...);
};
template <typename T>
using c_function_pointer_t = typename c_function_pointer<T>::type;


// BINDING
#define BIND(fn) [this]() { this->fn(); }
template <class C, typename R = void, typename I> _requires((!std::is_void_v<I>))
constexpr std::function<R(I)> bind(C* _this, R(C::* fn)(I)) {
	return [_this, fn](I data) -> R { return (_this->*fn)(data); };
}
template <class C, typename R = void>
constexpr std::function<R(void)> bind(C* _this, R(C::* fn)()) {
	return [_this, fn]() -> R { return (_this->*fn)(); };
}


// VISITOR
#ifdef HAS_CPP17
template <class... Ts>
struct visitor : Ts... { using Ts::operator()...; };
#endif


// VARIANT
#ifdef HAS_CPP17
// The requires are not necessary, but help with error messages, pointing directly to the function call, rather then in here
template <typename T, typename... VTs> _requires(is_one_of<T, VTs...>)
bool inline constexpr holds(const std::variant<VTs...>& var) { return std::holds_alternative<T>(var); }
template <typename T, typename... Vs> _requires((variant_contains<T, Vs>::value && ...))
bool inline constexpr holds(const Vs&... vars) { return (... && std::holds_alternative<T>(vars)); }
#endif


// CHRONO
void print_time(const stdc::nanoseconds& time, std::ostream& out = std::cout);


// FILE
#ifdef HAS_CPP17
std::string ReadFile(const fs::path& file);
#endif


struct Time
{
	Time(stdc::nanoseconds* res) {
		Res = res;
		Start = stdclock::now();
	}
	~Time() {
		*Res = stdclock::now() - Start;
	}

	stdc::nanoseconds* Res = nullptr;
	stdc::time_point<stdclock> Start;
};


struct AxisRange { float min, max; };
struct GraphRange { AxisRange x, y; };


// ################################################################## LOGGING ##################################################################
template <typename T> /*requires std::is_arithmetic_v<T>*/
std::ostream& operator<<(std::ostream& out, const std::vector<T>& data)
{
	out << '[';
	for (size_t i = 0; i + 1 < data.size(); i++)
		out << data[i] << ", ";
	if (!data.empty())
		out << data.back();
	out << ']';
	return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const coord<T, 2>& c) {
	return out << '[' << c.x << ':' << c.y << ']';
}
template <typename T>
std::ostream& operator<<(std::ostream& out, const coord<T, 3>& c) {
	return out << '[' << c.x << ':' << c.y << ':' << c.z << ']';
}
template <typename T>
std::ostream& operator<<(std::ostream& out, const coord<T, 4>& c) {
	return out << '[' << c.x << ':' << c.y << ':' << c.z << ':' << c.w << ']';
}
// ################################################################## LOGGING ##################################################################
