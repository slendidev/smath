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
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <smath/smath.hpp>

namespace smath
{

template<class Scalar, int N, int Options, int MaxRows, int MaxCols>
requires(std::is_arithmetic_v<Scalar> && (N > 0))
struct interop_adapter<Eigen::Matrix<Scalar, N, 1, Options, MaxRows, MaxCols>>
{
	using smath_type = Vec<static_cast<std::size_t>(N), Scalar>;

	static constexpr auto to_smath(
	    Eigen::Matrix<Scalar, N, 1, Options, MaxRows, MaxCols> const &v)
	    -> smath_type
	{
		smath_type out {};
		for (std::size_t i = 0; i < static_cast<std::size_t>(N); ++i)
			out[i] = v(static_cast<int>(i));
		return out;
	}

	static constexpr auto from_smath(smath_type const &v)
	    -> Eigen::Matrix<Scalar, N, 1, Options, MaxRows, MaxCols>
	{
		Eigen::Matrix<Scalar, N, 1, Options, MaxRows, MaxCols> out;
		for (std::size_t i = 0; i < static_cast<std::size_t>(N); ++i)
			out(static_cast<int>(i)) = v[i];
		return out;
	}
};

template<class Scalar, int N, int Options, int MaxRows, int MaxCols>
requires(std::is_arithmetic_v<Scalar> && (N > 0))
struct interop_adapter<Eigen::Matrix<Scalar, 1, N, Options, MaxRows, MaxCols>>
{
	using smath_type = Vec<static_cast<std::size_t>(N), Scalar>;

	static constexpr auto to_smath(
	    Eigen::Matrix<Scalar, 1, N, Options, MaxRows, MaxCols> const &v)
	    -> smath_type
	{
		smath_type out {};
		for (std::size_t i = 0; i < static_cast<std::size_t>(N); ++i)
			out[i] = v(static_cast<int>(i));
		return out;
	}

	static constexpr auto from_smath(smath_type const &v)
	    -> Eigen::Matrix<Scalar, 1, N, Options, MaxRows, MaxCols>
	{
		Eigen::Matrix<Scalar, 1, N, Options, MaxRows, MaxCols> out;
		for (std::size_t i = 0; i < static_cast<std::size_t>(N); ++i)
			out(static_cast<int>(i)) = v[i];
		return out;
	}
};

template<class Scalar, int R, int C, int Options, int MaxRows, int MaxCols>
requires(std::is_arithmetic_v<Scalar> && (R > 1) && (C > 1))
struct interop_adapter<Eigen::Matrix<Scalar, R, C, Options, MaxRows, MaxCols>>
{
	using smath_type
	    = Mat<static_cast<std::size_t>(R), static_cast<std::size_t>(C), Scalar>;

	static constexpr auto to_smath(
	    Eigen::Matrix<Scalar, R, C, Options, MaxRows, MaxCols> const &m)
	    -> smath_type
	{
		smath_type out {};
		for (std::size_t c = 0; c < static_cast<std::size_t>(C); ++c) {
			for (std::size_t r = 0; r < static_cast<std::size_t>(R); ++r)
				out[r, c] = m(static_cast<int>(r), static_cast<int>(c));
		}
		return out;
	}

	static constexpr auto from_smath(smath_type const &m)
	    -> Eigen::Matrix<Scalar, R, C, Options, MaxRows, MaxCols>
	{
		Eigen::Matrix<Scalar, R, C, Options, MaxRows, MaxCols> out;
		for (std::size_t c = 0; c < static_cast<std::size_t>(C); ++c) {
			for (std::size_t r = 0; r < static_cast<std::size_t>(R); ++r)
				out(static_cast<int>(r), static_cast<int>(c)) = m[r, c];
		}
		return out;
	}
};

template<class Scalar, int Options>
requires std::is_arithmetic_v<Scalar>
struct interop_adapter<Eigen::Quaternion<Scalar, Options>>
{
	using smath_type = Quaternion<Scalar>;

	static constexpr auto to_smath(Eigen::Quaternion<Scalar, Options> const &q)
	    -> smath_type
	{
		return { q.x(), q.y(), q.z(), q.w() };
	}

	static constexpr auto from_smath(smath_type const &q)
	    -> Eigen::Quaternion<Scalar, Options>
	{
		return { q.w(), q.x(), q.y(), q.z() };
	}
};

} // namespace smath
