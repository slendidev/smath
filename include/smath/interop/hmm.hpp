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

#include <HandmadeMath.h>

#include <smath/smath.hpp>

namespace smath
{

template<> struct interop_adapter<HMM_Vec2>
{
	using smath_type = Vec<2, float>;

	static constexpr auto to_smath(HMM_Vec2 const &v) -> smath_type
	{
		return { v.X, v.Y };
	}

	static constexpr auto from_smath(smath_type const &v) -> HMM_Vec2
	{
		return { v.x(), v.y() };
	}
};

template<> struct interop_adapter<HMM_Vec3>
{
	using smath_type = Vec<3, float>;

	static constexpr auto to_smath(HMM_Vec3 const &v) -> smath_type
	{
		return { v.X, v.Y, v.Z };
	}

	static constexpr auto from_smath(smath_type const &v) -> HMM_Vec3
	{
		return { v.x(), v.y(), v.z() };
	}
};

template<> struct interop_adapter<HMM_Vec4>
{
	using smath_type = Vec<4, float>;

	static constexpr auto to_smath(HMM_Vec4 const &v) -> smath_type
	{
		return { v.X, v.Y, v.Z, v.W };
	}

	static constexpr auto from_smath(smath_type const &v) -> HMM_Vec4
	{
		return { v.x(), v.y(), v.z(), v.w() };
	}
};

template<> struct interop_adapter<HMM_Mat4>
{
	using smath_type = Mat<4, 4, float>;

	static constexpr auto to_smath(HMM_Mat4 const &m) -> smath_type
	{
		smath_type out {};
		out[0, 0] = m.Columns[0].X;
		out[1, 0] = m.Columns[0].Y;
		out[2, 0] = m.Columns[0].Z;
		out[3, 0] = m.Columns[0].W;
		out[0, 1] = m.Columns[1].X;
		out[1, 1] = m.Columns[1].Y;
		out[2, 1] = m.Columns[1].Z;
		out[3, 1] = m.Columns[1].W;
		out[0, 2] = m.Columns[2].X;
		out[1, 2] = m.Columns[2].Y;
		out[2, 2] = m.Columns[2].Z;
		out[3, 2] = m.Columns[2].W;
		out[0, 3] = m.Columns[3].X;
		out[1, 3] = m.Columns[3].Y;
		out[2, 3] = m.Columns[3].Z;
		out[3, 3] = m.Columns[3].W;
		return out;
	}

	static constexpr auto from_smath(smath_type const &m) -> HMM_Mat4
	{
		HMM_Mat4 out {};
		out.Columns[0] = { m[0, 0], m[1, 0], m[2, 0], m[3, 0] };
		out.Columns[1] = { m[0, 1], m[1, 1], m[2, 1], m[3, 1] };
		out.Columns[2] = { m[0, 2], m[1, 2], m[2, 2], m[3, 2] };
		out.Columns[3] = { m[0, 3], m[1, 3], m[2, 3], m[3, 3] };
		return out;
	}
};

template<> struct interop_adapter<HMM_Quat>
{
	using smath_type = Quaternion<float>;

	static constexpr auto to_smath(HMM_Quat const &q) -> smath_type
	{
		return { q.X, q.Y, q.Z, q.W };
	}

	static constexpr auto from_smath(smath_type const &q) -> HMM_Quat
	{
		return { q.x(), q.y(), q.z(), q.w() };
	}
};

} // namespace smath
