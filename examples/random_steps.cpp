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

#include <print>
#include <random>

#include <smath.hpp>

auto main() -> int {
  using namespace smath;

  Vec2d point;
  std::random_device rd;
  std::mt19937 rng{rd()};
  std::uniform_real_distribution<> dis(-5, 5);
  int i = 0;
  do {
    Vec2d add{dis(rng), dis(rng)};
    auto const n = point + add;
    std::println("{}: {:.2f} + {:.2f} -> {:.2f}", i, point, add, n);
    point = n;

    i++;
  } while (i < 15);
}
