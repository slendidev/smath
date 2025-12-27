/*
 * smath - Single-file linear algebra math library for C++23.
 *
 * Copyright 2025 Slendi <slendi@socopon.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * You can define the following macros to change functionality:
 * - SMATH_IMPLICIT_CONVERSIONS
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <numbers>
#include <optional>
#include <type_traits>

#ifndef SMATH_ANGLE_UNIT
#	define SMATH_ANGLE_UNIT rad
#endif // SMATH_ANGLE_UNIT

namespace smath {

template<std::size_t N, typename T>
requires std::is_arithmetic_v<T> struct Vec;

namespace detail {

#define SMATH_STR(x) #x
#define SMATH_XSTR(x) SMATH_STR(x)

consteval bool streq(const char *a, const char *b)
{
	for (;; ++a, ++b) {
		if (*a != *b)
			return false;
		if (*a == '\0')
			return true;
	}
}

enum class AngularUnit {
	Radians,
	Degrees,
	Turns,
};

consteval std::optional<AngularUnit> parse_unit(char const *s)
{
	if (streq(s, "rad"))
		return AngularUnit::Radians;
	if (streq(s, "deg"))
		return AngularUnit::Degrees;
	if (streq(s, "turns"))
		return AngularUnit::Turns;
	return std::nullopt;
}

constexpr auto SMATH_ANGLE_UNIT_ID = parse_unit(SMATH_XSTR(SMATH_ANGLE_UNIT));
static_assert(SMATH_ANGLE_UNIT_ID != std::nullopt,
    "Invalid SMATH_ANGLE_UNIT. Should be rad, deg, or turns.");

template<std::size_t N> struct FixedString {
	char data[N] {};
	static constexpr std::size_t size = N - 1;
	constexpr FixedString(char const (&s)[N])
	{
		for (std::size_t i = 0; i < N; ++i)
			data[i] = s[i];
	}
	constexpr char operator[](std::size_t i) const { return data[i]; }
};
template<class X> struct is_Vec : std::false_type { };
template<std::size_t M, class U> struct is_Vec<Vec<M, U>> : std::true_type { };
template<class X>
inline constexpr bool is_Vec_v = is_Vec<std::remove_cvref_t<X>>::value;
template<class X>
inline constexpr bool is_scalar_v
    = std::is_arithmetic_v<std::remove_cvref_t<X>>;
template<class X> struct Vec_size;
template<std::size_t M, class U>
struct Vec_size<Vec<M, U>> : std::integral_constant<std::size_t, M> { };

template<class T> constexpr auto pack_unorm8(T v) -> std::uint8_t
{
	static_assert(std::is_floating_point_v<T>);
	T c = std::clamp(v, T(0), T(1));
	T scaled = c * T(255);
	int i = static_cast<int>(scaled + T(0.5));
	if (i < 0)
		i = 0;
	if (i > 255)
		i = 255;
	return static_cast<std::uint8_t>(i);
}

template<class T> constexpr auto pack_snorm8(T v) -> std::int8_t
{
	static_assert(std::is_floating_point_v<T>);
	T c = std::clamp(v, T(-1), T(1));
	T scaled = c * T(127);
	int i
	    = static_cast<int>(scaled >= T(0) ? scaled + T(0.5) : scaled - T(0.5));
	if (i < -127)
		i = -127;
	if (i > 127)
		i = 127;
	return static_cast<std::int8_t>(i);
}

template<class T> constexpr auto unpack_unorm8(std::uint8_t b) -> T
{
	static_assert(std::is_floating_point_v<T>);
	return static_cast<T>(b) / T(255);
}

template<class T> constexpr auto unpack_snorm8(std::int8_t b) -> T
{
	static_assert(std::is_floating_point_v<T>);
	int i = static_cast<int>(b);
	if (i < -127)
		i = -127;
	if (i > 127)
		i = 127;
	return static_cast<T>(i) / T(127);
}

} // namespace detail

template<std::size_t N, typename T = float>
requires std::is_arithmetic_v<T> struct Vec : std::array<T, N> {
private:
	template<class X> static consteval std::size_t extent()
	{
		if constexpr (detail::is_Vec_v<X>)
			return detail::Vec_size<std::remove_cvref_t<X>>::value;
		else if constexpr (detail::is_scalar_v<X>)
			return 1;
		else
			return 0; // Should be unreachable
	}
	template<class... Args> static consteval std::size_t total_extent()
	{
		return (extent<Args>() + ... + 0);
	}

public:
	// Constructors
	constexpr Vec() noexcept
	{
		for (auto &v : *this)
			v = T(0);
	}

	explicit constexpr Vec(T const &s) noexcept
	{
		for (auto &v : *this)
			v = s;
	}

	template<typename... Args>
	requires((detail::is_scalar_v<Args> || detail::is_Vec_v<Args>) && ...)
	    && (total_extent<Args...>() == N)
	    && (!(sizeof...(Args) == 1 && (detail::is_Vec_v<Args> && ...)))
	constexpr Vec(Args &&...args) noexcept
	{
		std::size_t i = 0;
		(fill_one(i, std::forward<Args>(args)), ...);
	}

	// Member accesses
	// NOTE: This can (probably) be improved with C++26 reflection in the
	// future.
#define VEC_ACC(component, req, idx) \
	constexpr auto component() noexcept -> T &requires(N >= req) { \
		return (*this)[idx]; \
	} constexpr auto component() const->T const & \
	requires(N >= req) \
	{ \
		return (*this)[idx]; \
	}

	VEC_ACC(r, 1, 0)
	VEC_ACC(g, 2, 1)
	VEC_ACC(b, 3, 2)
	VEC_ACC(a, 4, 3)

	VEC_ACC(x, 1, 0)
	VEC_ACC(y, 2, 1)
	VEC_ACC(z, 3, 2)
	VEC_ACC(w, 4, 3)

	VEC_ACC(s, 1, 0)
	VEC_ACC(t, 2, 1)
	VEC_ACC(p, 3, 2)
	VEC_ACC(q, 4, 3)

	VEC_ACC(u, 1, 0)
	VEC_ACC(v, 2, 1)
#undef VEC_ACC

	template<class... Args, std::size_t... Is>
	constexpr void unpack_impl(
	    std::index_sequence<Is...>, Args &...args) noexcept
	{
		((args = (*this)[Is]), ...);
	}

	template<class... Args> constexpr void unpack(Args &...args) noexcept
	{
		unpack_impl(std::index_sequence_for<Args...> {}, args...);
	}

	// Unary
	constexpr auto operator-() noexcept -> Vec
	{
		Vec r {};
		for (std::size_t i = 0; i < N; ++i)
			r[i] = -(*this)[i];
		return r;
	}

	// RHS operations
	friend constexpr auto operator+(T s, Vec const &v) noexcept -> Vec
	{
		return v + s;
	}
	friend constexpr auto operator-(T s, Vec const &v) noexcept -> Vec
	{
		return Vec(s) - v;
	}
	friend constexpr auto operator*(T s, Vec const &v) noexcept -> Vec
	{
		return v * s;
	}
	friend constexpr auto operator/(T s, Vec const &v) noexcept -> Vec
	{
		Vec r {};
		for (std::size_t i = 0; i < N; ++i)
			r[i] = s / v[i];
		return r;
	}

	// Members
#define VEC_OP(op) \
	constexpr auto operator op(Vec const &rhs) const noexcept -> Vec \
	{ \
		Vec result {}; \
		for (std::size_t i = 0; i < N; ++i) { \
			result[i] = (*this)[i] op rhs[i]; \
		} \
		return result; \
	} \
	constexpr auto operator op(T const &rhs) const noexcept -> Vec \
	{ \
		Vec result {}; \
		for (std::size_t i = 0; i < N; ++i) { \
			result[i] = (*this)[i] op rhs; \
		} \
		return result; \
	}
	VEC_OP(+)
	VEC_OP(-)
	VEC_OP(*)
	VEC_OP(/)
#undef VEC_OP
#define VEC_OP_ASSIGN(sym) \
	constexpr Vec &operator sym##=(Vec const &rhs) noexcept \
	{ \
		for (std::size_t i = 0; i < N; ++i) \
			(*this)[i] sym## = rhs[i]; \
		return *this; \
	} \
	constexpr Vec &operator sym##=(T const &s) noexcept \
	{ \
		for (std::size_t i = 0; i < N; ++i) \
			(*this)[i] sym## = s; \
		return *this; \
	}
	VEC_OP_ASSIGN(+)
	VEC_OP_ASSIGN(-)
	VEC_OP_ASSIGN(*)
	VEC_OP_ASSIGN(/)
#undef VEC_OP_ASSIGN

	constexpr auto operator==(Vec const &v) const noexcept -> bool
	{
		for (std::size_t i = 0; i < N; ++i)
			if ((*this)[i] != v[i])
				return false;
		return true;
	}

	constexpr auto operator!=(Vec const &v) const noexcept -> bool
	{
		return !(*this == v);
	}

	constexpr auto magnitude() const noexcept -> T
	requires std::is_floating_point_v<T>
	{
		T total = 0;
		for (auto const &v : *this)
			total += v * v;
		return std::sqrt(total);
	}
	constexpr auto length() const noexcept -> T
	requires std::is_floating_point_v<T>
	{
		return this->magnitude();
	}

	template<typename U = T>
	requires std::is_floating_point_v<U>
	constexpr auto normalized_safe(U eps = EPS_DEFAULT) const noexcept -> Vec
	{
		auto m = magnitude();
		return (m > eps) ? (*this) / m : Vec {};
	}
	template<typename U = T>
	requires std::is_floating_point_v<U>
	constexpr auto normalize_safe(U eps = EPS_DEFAULT) const noexcept -> Vec
	{
		return normalized_safe(eps);
	}

	[[nodiscard]] constexpr auto normalized() noexcept -> Vec<N, T>
	requires std::is_floating_point_v<T>
	{
		return (*this) / this->magnitude();
	}
	[[nodiscard]] constexpr auto normalize() noexcept -> Vec<N, T>
	requires std::is_floating_point_v<T>
	{
		return this->normalized();
	}
	[[nodiscard]] constexpr auto unit() noexcept -> Vec<N, T>
	requires std::is_floating_point_v<T>
	{
		return this->normalized();
	}

	[[nodiscard]] constexpr auto dot(Vec<N, T> const &other) const noexcept -> T
	{
		T res = 0;
		for (std::size_t i = 0; i < N; ++i) {
			res += (*this)[i] * other[i];
		}
		return res;
	}

	static constexpr T EPS_DEFAULT = T(1e-6);
	template<class U = T>
	requires std::is_floating_point_v<U>
	[[nodiscard]] constexpr auto approx_equal(
	    Vec const &rhs, U eps = EPS_DEFAULT) const noexcept
	{
		using F = std::conditional_t<std::is_floating_point_v<U>, U, double>;
		for (size_t i = 0; i < N; ++i)
			if (std::abs(F((*this)[i] - rhs[i])) > F(eps))
				return false;
		return true;
	}

	template<class U = T>
	constexpr auto magnitude_promoted() const noexcept
	requires std::is_floating_point_v<T>
	{
		using F = std::conditional_t<std::is_floating_point_v<U>, U, double>;
		F s = 0;
		for (auto v : *this)
			s += F(v) * F(v);
		return std::sqrt(s);
	}

	template<typename U = T>
	requires(N == 3) constexpr auto cross(Vec const &r) const noexcept -> Vec
	{
		return { (*this)[1] * r[2] - (*this)[2] * r[1],
			(*this)[2] * r[0] - (*this)[0] * r[2],
			(*this)[0] * r[1] - (*this)[1] * r[0] };
	}

	constexpr auto distance(Vec const &r) const noexcept -> T
	requires std::is_floating_point_v<T>
	{
		return (*this - r).magnitude();
	}

	constexpr auto project_onto(Vec const &n) const noexcept -> Vec
	requires std::is_floating_point_v<T>
	{
		auto d = this->dot(n);
		auto nn = n.dot(n);
		return (nn ? (d / nn) * n : Vec());
	}

	template<class U>
	requires(std::is_arithmetic_v<U> && N >= 1)
	constexpr explicit(!std::is_convertible_v<U, T>)
	    Vec(Vec<N, U> const &other) noexcept
	{
		for (std::size_t i = 0; i < N; ++i)
			this->operator[](i) = static_cast<T>(other[i]);
	}

	template<class U>
	requires(std::is_arithmetic_v<U> && N >= 1) constexpr explicit(
	    !std::is_convertible_v<T, U>) operator Vec<N, U>() const noexcept
	{
		Vec<N, U> r {};
		for (std::size_t i = 0; i < N; ++i)
			r[i] = static_cast<U>((*this)[i]);
		return r;
	}

	template<class U>
	requires(std::is_arithmetic_v<U> && !std::is_same_v<U, T>)
	constexpr auto operator=(Vec<N, U> const &rhs) noexcept -> Vec &
	{
		for (std::size_t i = 0; i < N; ++i)
			(*this)[i] = static_cast<T>(rhs[i]);
		return *this;
	}

private:
	constexpr void fill_one(std::size_t &i, T const &v) noexcept
	{
		(*this)[i++] = v;
	}
#ifdef SMATH_IMPLICIT_CONVERSIONS
	template<class U>
	requires std::is_arithmetic_v<U>
	    && (!std::is_same_v<U, T>)constexpr void fill_one(
	    std::size_t &i, const U &v) noexcept
	{
		(*this)[i++] = static_cast<T>(v);
	}
	template<std::size_t M, class U>
	constexpr void fill_one(std::size_t &i, Vec<M, U> const &v) noexcept
	{
		for (std::size_t k = 0; k < M; ++k)
			(*this)[i++] = static_cast<T>(v[k]);
	}
#endif // SMATH_IMPLICIT_CONVERSIONS
	template<std::size_t M>
	constexpr void fill_one(std::size_t &i, const Vec<M, T> &v) noexcept
	{
		for (std::size_t k = 0; k < M; ++k)
			(*this)[i++] = static_cast<T>(v[k]);
	}
};

template<size_t I, size_t N, class T> constexpr T &get(Vec<N, T> &v) noexcept
{
	static_assert(I < N);
	return v[I];
}
template<size_t I, size_t N, class T>
constexpr T const &get(Vec<N, T> const &v) noexcept
{
	static_assert(I < N);
	return v[I];
}
template<size_t I, size_t N, class T> constexpr T &&get(Vec<N, T> &&v) noexcept
{
	static_assert(I < N);
	return std::move(v[I]);
}
template<size_t I, size_t N, class T>
constexpr T const &&get(Vec<N, T> const &&v) noexcept
{
	static_assert(I < N);
	return std::move(v[I]);
}

template<std::size_t N, typename T = float>
requires std::is_arithmetic_v<T>
using VecOrScalar = std::conditional_t<N == 1, T, Vec<N, T>>;

namespace detail {

consteval auto char_to_idx(char c) -> std::size_t
{
	if (c == 'r' || c == 'x' || c == 's' || c == 'u')
		return 0;
	else if (c == 'g' || c == 'y' || c == 't' || c == 'v')
		return 1;
	else if (c == 'b' || c == 'z' || c == 'p')
		return 2;
	else if (c == 'a' || c == 'w' || c == 'q')
		return 3;
	return static_cast<std::size_t>(-1);
}

constexpr auto is_valid(char c) -> bool
{
	switch (c) {
	case 'r':
	case 'g':
	case 'b':
	case 'a':
	case 'x':
	case 'y':
	case 'z':
	case 'w':
	case 's':
	case 't':
	case 'p':
	case 'q':
	case 'u':
	case 'v':
		return true;
	default:
		return false;
	}
	return false;
}

template<detail::FixedString S, std::size_t N, typename T, std::size_t... I>
constexpr auto swizzle_impl(Vec<N, T> const &v, std::index_sequence<I...>)
    -> VecOrScalar<S.size, T>
{
	static_assert(((is_valid(S[I])) && ...), "Invalid swizzle component");
	static_assert(
	    ((char_to_idx(S[I]) < N) && ...), "Pattern index out of bounds");
	VecOrScalar<S.size, T> out {};
	std::size_t i = 0;
	((out[i++] = v[char_to_idx(S[I])]), ...);
	return out;
}

template<FixedString S>
concept SwizzleCharsOK = []<std::size_t... I>(std::index_sequence<I...>) {
	return ((is_valid(S[I])) && ...);
}(std::make_index_sequence<S.size> {});

template<FixedString S, std::size_t N>
concept SwizzleInBounds = []<std::size_t... I>(std::index_sequence<I...>) {
	return ((char_to_idx(S[I]) < N) && ...);
}(std::make_index_sequence<S.size> {});

template<FixedString S, std::size_t N>
concept ValidSwizzle
    = (S.size > 0) && SwizzleCharsOK<S> && SwizzleInBounds<S, N>;

} // namespace detail

template<detail::FixedString S, std::size_t N, typename T>
requires detail::ValidSwizzle<S, N>
constexpr auto swizzle(Vec<N, T> const &v) -> VecOrScalar<S.size, T>
{
	return detail::swizzle_impl<S>(v, std::make_index_sequence<S.size> {});
}

using Vec2 = Vec<2>;
using Vec3 = Vec<3>;
using Vec4 = Vec<4>;

using Vec2d = Vec<2, double>;
using Vec3d = Vec<3, double>;
using Vec4d = Vec<4, double>;

template<class T> constexpr auto deg(T const value)
{
	using R = std::common_type_t<T, double>;
	if constexpr (detail::SMATH_ANGLE_UNIT_ID == detail::AngularUnit::Degrees) {
		return static_cast<R>(value);
	} else if constexpr (detail::SMATH_ANGLE_UNIT_ID
	    == detail::AngularUnit::Radians) {
		return static_cast<R>(value) * static_cast<R>(std::numbers::pi / 180.0);
	} else if constexpr (detail::SMATH_ANGLE_UNIT_ID
	    == detail::AngularUnit::Turns) {
		return static_cast<R>(value) / static_cast<R>(360.0);
	}
}

template<class T> constexpr auto rad(T const value)
{
	using R = std::common_type_t<T, double>;
	if constexpr (detail::SMATH_ANGLE_UNIT_ID == detail::AngularUnit::Degrees) {
		return static_cast<R>(value) * static_cast<R>(180.0 / std::numbers::pi);
	} else if constexpr (detail::SMATH_ANGLE_UNIT_ID
	    == detail::AngularUnit::Radians) {
		return static_cast<R>(value);
	} else if constexpr (detail::SMATH_ANGLE_UNIT_ID
	    == detail::AngularUnit::Turns) {
		return static_cast<R>(value)
		    / (static_cast<R>(2.0) * static_cast<R>(std::numbers::pi));
	}
}

template<class T> constexpr auto turns(T const value)
{
	using R = std::common_type_t<T, double>;
	if constexpr (detail::SMATH_ANGLE_UNIT_ID == detail::AngularUnit::Degrees) {
		return static_cast<R>(value) * static_cast<R>(360.0);
	} else if constexpr (detail::SMATH_ANGLE_UNIT_ID
	    == detail::AngularUnit::Radians) {
		return static_cast<R>(value)
		    * (static_cast<R>(2.0) * static_cast<R>(std::numbers::pi));
	} else if constexpr (detail::SMATH_ANGLE_UNIT_ID
	    == detail::AngularUnit::Turns) {
		return static_cast<R>(value);
	}
}

template<class T> struct Quaternion : Vec<4, T> {
	using Base = Vec<4, T>;
	using Base::Base;
	using Base::operator=;

	constexpr Base &vec() noexcept { return *this; }
	constexpr Base const &vec() const noexcept { return *this; }

	constexpr T &x() noexcept { return Base::x(); }
	constexpr T const &x() const noexcept { return Base::x(); }
	constexpr T &y() noexcept { return Base::y(); }
	constexpr T const &y() const noexcept { return Base::y(); }
	constexpr T &z() noexcept { return Base::z(); }
	constexpr T const &z() const noexcept { return Base::z(); }
	constexpr T &w() noexcept { return Base::w(); }
	constexpr T const &w() const noexcept { return Base::w(); }

	constexpr auto operator*(Quaternion const &rhs) const noexcept -> Quaternion
	{
		Quaternion r;
		auto const &a = *this;

		r.x() = a.w() * rhs.x() + a.x() * rhs.w() + a.y() * rhs.z()
		    - a.z() * rhs.y();
		r.y() = a.w() * rhs.y() - a.x() * rhs.z() + a.y() * rhs.w()
		    + a.z() * rhs.x();
		r.z() = a.w() * rhs.z() + a.x() * rhs.y() - a.y() * rhs.x()
		    + a.z() * rhs.w();
		r.w() = a.w() * rhs.w() - a.x() * rhs.x() - a.y() * rhs.y()
		    - a.z() * rhs.z();

		return r;
	}
};

template<class T>
requires std::is_floating_point_v<T>
[[nodiscard]] constexpr auto pack_unorm4x8(Vec<4, T> const &v) -> std::uint32_t
{
	std::uint32_t r = detail::pack_unorm8(v[0]);
	std::uint32_t g = detail::pack_unorm8(v[1]);
	std::uint32_t b = detail::pack_unorm8(v[2]);
	std::uint32_t a = detail::pack_unorm8(v[3]);

	return (r) | (g << 8) | (b << 16) | (a << 24);
}

template<class T>
requires std::is_floating_point_v<T>
[[nodiscard]] constexpr auto pack_snorm4x8(Vec<4, T> const &v) -> std::uint32_t
{
	std::uint32_t r = static_cast<std::uint8_t>(detail::pack_snorm8(v[0]));
	std::uint32_t g = static_cast<std::uint8_t>(detail::pack_snorm8(v[1]));
	std::uint32_t b = static_cast<std::uint8_t>(detail::pack_snorm8(v[2]));
	std::uint32_t a = static_cast<std::uint8_t>(detail::pack_snorm8(v[3]));

	return (r) | (g << 8) | (b << 16) | (a << 24);
}

template<class T = float>
requires std::is_floating_point_v<T>
[[nodiscard]] constexpr auto unpack_unorm4x8(std::uint32_t packed) -> Vec<4, T>
{
	std::uint8_t r = static_cast<std::uint8_t>(packed & 0xFFu);
	std::uint8_t g = static_cast<std::uint8_t>((packed >> 8) & 0xFFu);
	std::uint8_t b = static_cast<std::uint8_t>((packed >> 16) & 0xFFu);
	std::uint8_t a = static_cast<std::uint8_t>((packed >> 24) & 0xFFu);

	return {
		detail::unpack_unorm8<T>(r),
		detail::unpack_unorm8<T>(g),
		detail::unpack_unorm8<T>(b),
		detail::unpack_unorm8<T>(a),
	};
}

template<class T = float>
requires std::is_floating_point_v<T>
[[nodiscard]] constexpr auto unpack_snorm4x8(std::uint32_t packed) -> Vec<4, T>
{
	std::int8_t r = static_cast<std::int8_t>(packed & 0xFFu);
	std::int8_t g = static_cast<std::int8_t>((packed >> 8) & 0xFFu);
	std::int8_t b = static_cast<std::int8_t>((packed >> 16) & 0xFFu);
	std::int8_t a = static_cast<std::int8_t>((packed >> 24) & 0xFFu);

	return {
		detail::unpack_snorm8<T>(r),
		detail::unpack_snorm8<T>(g),
		detail::unpack_snorm8<T>(b),
		detail::unpack_snorm8<T>(a),
	};
}

template<std::size_t R, std::size_t C, typename T = float>
requires std::is_arithmetic_v<T> struct Mat : std::array<Vec<R, T>, C> {
	using Base = std::array<Vec<R, T>, C>;
	using Base::operator[];

	constexpr auto operator[](std::size_t const row, std::size_t const column)
	    -> T &
	{
		return col(column)[row];
	}
	constexpr auto operator[](
	    std::size_t const row, std::size_t const column) const -> T const &
	{
		return col(column)[row];
	}

	constexpr Mat() noexcept
	{
		for (auto &col : *this)
			col = Vec<R, T> {};
	}

	constexpr explicit Mat(T const &diag) noexcept
	requires(R == C)
	{
		for (std::size_t c = 0; c < C; ++c) {
			(*this)[c] = Vec<R, T> {};
			(*this)[c][c] = diag;
		}
	}

	template<typename... Cols>
	requires(sizeof...(Cols) == C
	    && (std::same_as<std::remove_cvref_t<Cols>, Vec<R, T>> && ...))
	constexpr Mat(Cols const &...cols) noexcept
	    : Base { cols... }
	{
	}

	constexpr auto col(std::size_t j) noexcept -> Vec<R, T> &
	{
		return (*this)[j];
	}
	constexpr auto col(std::size_t j) const noexcept -> Vec<R, T> const &
	{
		return (*this)[j];
	}

	constexpr auto operator()(std::size_t row, std::size_t col) noexcept -> T &
	{
		return (*this)[col][row];
	}
	constexpr auto operator()(std::size_t row, std::size_t col) const noexcept
	    -> T const &
	{
		return (*this)[col][row];
	}

	constexpr auto operator-() const noexcept -> Mat
	{
		Mat r {};
		for (std::size_t c = 0; c < C; ++c)
			r[c] = -(*this)[c];
		return r;
	}

	constexpr auto operator+=(Mat const &rhs) noexcept -> Mat &
	{
		for (std::size_t c = 0; c < C; ++c)
			(*this)[c] += rhs[c];
		return *this;
	}
	constexpr auto operator-=(Mat const &rhs) noexcept -> Mat &
	{
		for (std::size_t c = 0; c < C; ++c)
			(*this)[c] -= rhs[c];
		return *this;
	}
	friend constexpr auto operator+(Mat lhs, Mat const &rhs) noexcept -> Mat
	{
		lhs += rhs;
		return lhs;
	}
	friend constexpr auto operator-(Mat lhs, Mat const &rhs) noexcept -> Mat
	{
		lhs -= rhs;
		return lhs;
	}

	constexpr auto operator*=(T const &s) noexcept -> Mat &
	{
		for (std::size_t c = 0; c < C; ++c)
			(*this)[c] *= s;
		return *this;
	}
	constexpr auto operator/=(T const &s) noexcept -> Mat &
	{
		for (std::size_t c = 0; c < C; ++c)
			(*this)[c] /= s;
		return *this;
	}
	friend constexpr auto operator*(Mat lhs, T const &s) noexcept -> Mat
	{
		lhs *= s;
		return lhs;
	}
	friend constexpr auto operator*(T const &s, Mat rhs) noexcept -> Mat
	{
		rhs *= s;
		return rhs;
	}
	friend constexpr auto operator/(Mat lhs, T const &s) noexcept -> Mat
	{
		lhs /= s;
		return lhs;
	}

	[[nodiscard]] constexpr auto operator==(Mat const &rhs) const noexcept
	    -> bool
	{
		for (std::size_t c = 0; c < C; ++c)
			if (!((*this)[c] == rhs[c]))
				return false;
		return true;
	}
	[[nodiscard]] constexpr auto operator!=(Mat const &rhs) const noexcept
	    -> bool
	{
		return !(*this == rhs);
	}

	static constexpr T EPS_DEFAULT = T(1e-6);
	template<class U = T>
	requires std::is_floating_point_v<U>
	[[nodiscard]] constexpr auto approx_equal(
	    Mat const &rhs, U eps = EPS_DEFAULT) const noexcept -> bool
	{
		for (std::size_t c = 0; c < C; ++c)
			if (!(*this)[c].approx_equal(rhs[c], eps))
				return false;
		return true;
	}

	[[nodiscard]] constexpr auto transposed() const noexcept -> Mat<C, R, T>
	{
		Mat<C, R, T> r {};
		for (std::size_t c = 0; c < C; ++c)
			for (std::size_t r_idx = 0; r_idx < R; ++r_idx)
				r(r_idx, c) = (*this)(c, r_idx);
		return r;
	}

	[[nodiscard]] static constexpr auto identity() noexcept -> Mat<R, C, T>
	requires(R == C)
	{
		Mat<R, C, T> m {};
		for (std::size_t i = 0; i < R; ++i)
			m(i, i) = T(1);
		return m;
	}
};

using Mat2 = Mat<2, 2>;
using Mat3 = Mat<3, 3>;
using Mat4 = Mat<4, 4>;

using Mat2d = Mat<2, 2, double>;
using Mat3d = Mat<3, 3, double>;
using Mat4d = Mat<4, 4, double>;

template<std::size_t R, std::size_t C, typename T>
[[nodiscard]] constexpr Vec<R, T> operator*(
    Mat<R, C, T> const &m, Vec<C, T> const &v) noexcept
{
	Vec<R, T> out {};
	for (std::size_t c = 0; c < C; ++c)
		out += m.col(c) * v[c];
	return out;
}

// Matrix * Matrix
template<std::size_t R, std::size_t C, std::size_t K, typename T>
[[nodiscard]] constexpr Mat<R, K, T> operator*(
    Mat<R, C, T> const &a, Mat<C, K, T> const &b) noexcept
{
	Mat<R, K, T> out {};
	for (std::size_t k = 0; k < K; ++k) {
		for (std::size_t r = 0; r < R; ++r) {
			T sum = T(0);
			for (std::size_t c = 0; c < C; ++c)
				sum += a(r, c) * b(c, k);
			out(r, k) = sum;
		}
	}
	return out;
}

// Mat3 transformations
template<typename T>
[[nodiscard]] inline auto translate(Mat<3, 3, T> const &m, Vec<2, T> const &v)
    -> Mat<3, 3, T>
{
	Mat<3, 3, T> res { m };
	res[2] = m[0] * v[0] + m[1] * v[1] + m[2];
	return res;
}

template<typename T>
[[nodiscard]] inline auto translate(Vec<2, T> const &v) -> Mat<3, 3, T>
{
	Mat<3, 3, T> res { 1 };
	res[2].x() = v.x();
	res[2].y() = v.y();
	return res;
}

template<typename T>
[[nodiscard]] inline auto rotate(Mat<3, 3, T> const &m, T const angle)
    -> Mat<3, 3, T>
{
	Mat<3, 3, T> res;

	T const c { std::cos(angle) };
	T const s { std::sin(angle) };

	res[0] = m[0] * c + m[1] * s;
	res[1] = m[0] * -s + m[1] * c;
	res[2] = m[2];

	return res;
}

template<typename T>
[[nodiscard]] inline auto scale(Mat<3, 3, T> const &m, Vec<2, T> const &v)
    -> Mat<3, 3, T>
{
	Mat<3, 3, T> res;

	res[0] = m[0] * v[0];
	res[1] = m[1] * v[1];
	res[2] = m[2];

	return res;
}

template<typename T>
[[nodiscard]] inline auto shear_x(Mat<3, 3, T> const &m, T const v)
    -> Mat<3, 3, T>
{
	Mat<3, 3, T> res { 1 };
	res[1][0] = v;
	return m * res;
}

template<typename T>
[[nodiscard]] inline auto shear_y(Mat<3, 3, T> const &m, T const v)
    -> Mat<3, 3, T>
{
	Mat<3, 3, T> res { 1 };
	res[0][1] = v;
	return m * res;
}

// Mat4 transformations
template<typename T>
[[nodiscard]] inline auto translate(Mat<4, 4, T> const &m, Vec<3, T> const &v)
    -> Mat<4, 4, T>
{
	Mat<4, 4, T> res { m };
	res[3] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2] + m[3];
	return res;
}

template<typename T>
[[nodiscard]] inline auto translate(Vec<3, T> const &v) -> Mat<4, 4, T>
{
	Mat<4, 4, T> res { 1 };
	res[3].x() = v.x();
	res[3].y() = v.y();
	res[3].z() = v.z();
	return res;
}

template<typename T>
[[nodiscard]] inline auto rotate(Mat<4, 4, T> const &m, T const angle)
    -> Mat<4, 4, T>
{
	Mat<4, 4, T> res;

	T const c { std::cos(angle) };
	T const s { std::sin(angle) };

	res[0] = m[0] * c + m[1] * s;
	res[1] = m[0] * -s + m[1] * c;
	res[2] = m[2];
	res[3] = m[3];

	return res;
}

template<typename T>
[[nodiscard]] inline auto scale(Mat<4, 4, T> const &m, Vec<3, T> const &v)
    -> Mat<4, 4, T>
{
	Mat<4, 4, T> res;

	res[0] = m[0] * v[0];
	res[1] = m[1] * v[1];
	res[2] = m[2] * v[2];
	res[3] = m[3];

	return res;
}

template<typename T>
[[nodiscard]] inline auto shear_x(Mat<4, 4, T> const &m, T const v)
    -> Mat<4, 4, T>
{
	Mat<4, 4, T> res { 1 };
	res[0, 1] = v;
	return m * res;
}

template<typename T>
[[nodiscard]] inline auto shear_y(Mat<4, 4, T> const &m, T const v)
    -> Mat<4, 4, T>
{
	Mat<4, 4, T> res { 1 };
	res[1, 0] = v;
	return m * res;
}

template<typename T>
[[nodiscard]] inline auto shear_z(Mat<4, 4, T> const &m, T const v)
    -> Mat<4, 4, T>
{
	Mat<4, 4, T> res { 1 };
	res[2, 0] = v;
	return m * res;
}

template<typename T>
[[nodiscard]] inline auto matrix_ortho3d(T const left, T const right,
    T const bottom, T const top, T const near, T const far,
    bool const flip_z_axis = true) -> Mat<4, 4, T>
{
	Mat<4, 4, T> res {};

	res[0, 0] = 2 / (right - left);
	res[1, 1] = 2 / (top - bottom);
	res[2, 2] = -2 / (far - near);
	res[0, 3] = -(right + left) / (right - left);
	res[1, 3] = -(top + bottom) / (top - bottom);
	res[2, 3] = -(far + near) / (far - near);
	res[3, 3] = 1;

	if (flip_z_axis) {
		res[2] = -res[2];
	}

	return res;
}

template<typename T>
inline auto matrix_perspective(
    T fovy, T aspect, T znear, T zfar, bool flip_z_axis = false) -> Mat<4, 4, T>
{
	Mat<4, 4, T> m {};

	T const f { T(1) / std::tan(fovy / T(2)) };

	m[0, 0] = f / aspect;
	m[1, 1] = f;

	if (!flip_z_axis) {
		m[2, 2] = -(zfar + znear) / (zfar - znear);
		m[2, 3] = -(T(2) * zfar * znear) / (zfar - znear);
		m[3, 2] = -1;
		m[3, 3] = 0;
	} else {
		m[2, 2] = (zfar + znear) / (zfar - znear);
		m[2, 3] = (T(2) * zfar * znear) / (zfar - znear);
		m[3, 2] = 1;
		m[3, 3] = 0;
	}

	return m;
}

template<typename T>
[[nodiscard]] inline auto matrix_look_at(Vec<3, T> const eye,
    Vec<3, T> const center, Vec<3, T> const up, bool flip_z_axis = false)
    -> Mat<4, 4, T>
{
	auto f = (center - eye).normalized();
	auto s = f.cross(up).normalized();
	auto u = s.cross(f);

	if (!flip_z_axis) {
		return Mat<4, 4, T> {
			Vec<4, T> { s.x(), s.y(), s.z(), 0 },
			Vec<4, T> { u.x(), u.y(), u.z(), 0 },
			Vec<4, T> { -f.x(), -f.y(), -f.z(), 0 },
			Vec<4, T> { -s.dot(eye), -u.dot(eye), f.dot(eye), 1 },
		};
	} else {
		return Mat<4, 4, T> {
			Vec<4, T> { s.x(), s.y(), s.z(), 0 },
			Vec<4, T> { u.x(), u.y(), u.z(), 0 },
			Vec<4, T> { f.x(), f.y(), f.z(), 0 },
			Vec<4, T> { -s.dot(eye), -u.dot(eye), -f.dot(eye), 1 },
		};
	}
}

template<typename T>
[[nodiscard]] inline auto matrix_infinite_perspective(T const fovy,
    T const aspect, T const znear, bool flip_z_axis = false) -> Mat<4, 4, T>
{
	Mat<4, 4, T> m {};

	T const f = 1 / std::tan(fovy / T(2));
	m(0, 0) = f / aspect;
	m(1, 1) = f;

	if (!flip_z_axis) {
		m(2, 2) = -1;
		m(2, 3) = -T(2) * znear;
		m(3, 2) = -1;
	} else {
		m(2, 2) = 1;
		m(2, 3) = T(2) * znear;
		m(3, 2) = 1;
	}

	return m;
}

} // namespace smath

template<std::size_t N, typename T>
requires std::formattable<T, char>
struct std::formatter<smath::Vec<N, T>> : std::formatter<T> {
	constexpr auto parse(std::format_parse_context &ctx)
	{
		return std::formatter<T>::parse(ctx);
	}

	template<typename Ctx>
	auto format(smath::Vec<N, T> const &v, Ctx &ctx) const
	{
		auto out = ctx.out();
		*out++ = '{';
		for (std::size_t i = 0; i < N; ++i) {
			if (i) {
				*out++ = ',';
				*out++ = ' ';
			}
			out = std::formatter<T>::format(v[i], ctx);
		}
		*out++ = '}';
		return out;
	}
};

namespace std {
template<size_t N, class T>
struct tuple_size<smath::Vec<N, T>> : std::integral_constant<size_t, N> { };

template<size_t I, size_t N, class T>
struct tuple_element<I, smath::Vec<N, T>> {
	static_assert(I < N);
	using type = T;
};

} // namespace std
