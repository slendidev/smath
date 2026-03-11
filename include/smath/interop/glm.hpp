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

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include <smath/smath.hpp>

namespace smath
{

template<glm::length_t L, class T, glm::qualifier Q>
requires std::is_arithmetic_v<T>
struct interop_adapter<glm::vec<L, T, Q>>
{
	using smath_type = Vec<static_cast<std::size_t>(L), T>;

	static constexpr auto to_smath(glm::vec<L, T, Q> const &v) -> smath_type
	{
		smath_type out {};
		for (std::size_t i = 0; i < static_cast<std::size_t>(L); ++i)
			out[i] = v[static_cast<glm::length_t>(i)];
		return out;
	}

	static constexpr auto from_smath(smath_type const &v) -> glm::vec<L, T, Q>
	{
		glm::vec<L, T, Q> out {};
		for (std::size_t i = 0; i < static_cast<std::size_t>(L); ++i)
			out[static_cast<glm::length_t>(i)] = v[i];
		return out;
	}
};

template<glm::length_t C, glm::length_t R, class T, glm::qualifier Q>
requires std::is_arithmetic_v<T>
struct interop_adapter<glm::mat<C, R, T, Q>>
{
	using smath_type
	    = Mat<static_cast<std::size_t>(R), static_cast<std::size_t>(C), T>;

	static constexpr auto to_smath(glm::mat<C, R, T, Q> const &m) -> smath_type
	{
		smath_type out {};
		for (std::size_t c = 0; c < static_cast<std::size_t>(C); ++c) {
			for (std::size_t r = 0; r < static_cast<std::size_t>(R); ++r) {
				out[r, c] = m[static_cast<glm::length_t>(c)]
				             [static_cast<glm::length_t>(r)];
			}
		}
		return out;
	}

	static constexpr auto from_smath(smath_type const &m)
	    -> glm::mat<C, R, T, Q>
	{
		glm::mat<C, R, T, Q> out {};
		for (std::size_t c = 0; c < static_cast<std::size_t>(C); ++c) {
			for (std::size_t r = 0; r < static_cast<std::size_t>(R); ++r) {
				out[static_cast<glm::length_t>(c)]
				   [static_cast<glm::length_t>(r)]
				    = m[r, c];
			}
		}
		return out;
	}
};

template<class T, glm::qualifier Q>
requires std::is_arithmetic_v<T>
struct interop_adapter<glm::qua<T, Q>>
{
	using smath_type = Quaternion<T>;

	static constexpr auto to_smath(glm::qua<T, Q> const &q) -> smath_type
	{
		return { q.x, q.y, q.z, q.w };
	}

	static constexpr auto from_smath(smath_type const &q) -> glm::qua<T, Q>
	{
		return { q.w(), q.x(), q.y(), q.z() };
	}
};

} // namespace smath
