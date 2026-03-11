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

#include <SFML/System/Vector2.hpp>
#include <SFML/System/Vector3.hpp>

#include <smath/smath.hpp>

namespace smath
{

template<class T>
requires std::is_arithmetic_v<T>
struct interop_adapter<sf::Vector2<T>>
{
	using smath_type = Vec<2, T>;

	static constexpr auto to_smath(sf::Vector2<T> const &v) -> smath_type
	{
		return { v.x, v.y };
	}

	static constexpr auto from_smath(smath_type const &v) -> sf::Vector2<T>
	{
		return { v.x(), v.y() };
	}
};

template<class T>
requires std::is_arithmetic_v<T>
struct interop_adapter<sf::Vector3<T>>
{
	using smath_type = Vec<3, T>;

	static constexpr auto to_smath(sf::Vector3<T> const &v) -> smath_type
	{
		return { v.x, v.y, v.z };
	}

	static constexpr auto from_smath(smath_type const &v) -> sf::Vector3<T>
	{
		return { v.x(), v.y(), v.z() };
	}
};

} // namespace smath
