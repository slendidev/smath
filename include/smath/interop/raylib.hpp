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

#include <raylib.h>

#include <smath/smath.hpp>

namespace smath
{

template<> struct interop_adapter<Vector2>
{
	using smath_type = Vec<2, float>;

	static constexpr auto to_smath(Vector2 const &v) -> smath_type
	{
		return { v.x, v.y };
	}

	static constexpr auto from_smath(smath_type const &v) -> Vector2
	{
		return { v.x(), v.y() };
	}
};

template<> struct interop_adapter<Vector3>
{
	using smath_type = Vec<3, float>;

	static constexpr auto to_smath(Vector3 const &v) -> smath_type
	{
		return { v.x, v.y, v.z };
	}

	static constexpr auto from_smath(smath_type const &v) -> Vector3
	{
		return { v.x(), v.y(), v.z() };
	}
};

template<> struct interop_adapter<Vector4>
{
	using smath_type = Vec<4, float>;

	static constexpr auto to_smath(Vector4 const &v) -> smath_type
	{
		return { v.x, v.y, v.z, v.w };
	}

	static constexpr auto from_smath(smath_type const &v) -> Vector4
	{
		return { v.x(), v.y(), v.z(), v.w() };
	}
};

template<> struct interop_adapter<Matrix>
{
	using smath_type = Mat<4, 4, float>;

	static constexpr auto to_smath(Matrix const &m) -> smath_type
	{
		smath_type out {};
		out[0, 0] = m.m0;
		out[1, 0] = m.m1;
		out[2, 0] = m.m2;
		out[3, 0] = m.m3;
		out[0, 1] = m.m4;
		out[1, 1] = m.m5;
		out[2, 1] = m.m6;
		out[3, 1] = m.m7;
		out[0, 2] = m.m8;
		out[1, 2] = m.m9;
		out[2, 2] = m.m10;
		out[3, 2] = m.m11;
		out[0, 3] = m.m12;
		out[1, 3] = m.m13;
		out[2, 3] = m.m14;
		out[3, 3] = m.m15;
		return out;
	}

	static constexpr auto from_smath(smath_type const &m) -> Matrix
	{
		return {
			m[0, 0],
			m[0, 1],
			m[0, 2],
			m[0, 3],
			m[1, 0],
			m[1, 1],
			m[1, 2],
			m[1, 3],
			m[2, 0],
			m[2, 1],
			m[2, 2],
			m[2, 3],
			m[3, 0],
			m[3, 1],
			m[3, 2],
			m[3, 3],
		};
	}
};

} // namespace smath
