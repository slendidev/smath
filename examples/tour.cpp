/* Copyright 2025 Slendi <slendi@socopon.com>
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

#include <print>

// #define SMATH_IMPLICIT_CONVERSIONS

#include <smath.hpp>

int main() {
  using namespace smath;
  Vec3 v{1, 2, 3};
  std::println("v: {}", v);
  auto v2 = swizzle<"zyx">(v);
  std::println("v2: {}", v2);
  std::println("+: {}", v + v2);
  std::println("-: {}", v - v2);
  std::println("*: {}", v * v2);
  std::println("/: {}", v / v2);
  std::println("dot: {}", v.dot(v2));
  std::println("rrggbb: {}", swizzle<"rrggbb">(v));
  std::println("Magnitude: {}", v.magnitude());
  std::println("Normalized: {}", v.normalized());
  std::println("(alias) Unit: {}", v.unit());
  std::println("(alias) Normalize: {}", v.normalize());
  std::println("(alias) Length: {}", v.length());
  std::println("std::get<1>(v): {}", std::get<1>(v));
  auto [x, y, z] = v;
  std::println("Bindings: [{}, {}, {}]", x, y, z);
  float x1{}, y1{}, z1{};
  v.unpack(x1, y1, z1);
  std::println("Unpacked: {}, {}, {}", x1, y1, z1);

  // Let's mix and match!
  Vec<6> v3(v, 7, swizzle<"zy">(v2));
  std::println("{{v, 7, XZ(v2)}}: {}", v3);

  // Scalar operations
  std::println("v + 3: {}", v + 3);
  std::println("v - 3: {}", v - 3);
  std::println("v * 3: {}", v * 3);
  std::println("v / 3: {}", v / 3);

  std::println("3 + v: {}", 3 + v);
  std::println("3 - v: {}", 3 - v);
  std::println("3 * v: {}", 3 * v);
  std::println("3 / v: {}", 3 / v);

  // Casting
  auto v4 = static_cast<Vec3d>(v);
#ifdef SMATH_IMPLICIT_CONVERSIONS
  Vec3d v5 = v;
#else
  Vec3d v5 = static_cast<Vec3d>(v);
#endif // SMATH_IMPLICIT_CONVERSIONS
  std::println("Are v4 and v5 same? {}", v4 == v5);
}
