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

#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <type_traits>

namespace smath {

template <std::size_t N, typename T>
  requires std::is_arithmetic_v<T>
struct Vec;

namespace detail {

template <std::size_t N> struct FixedString {
  char data[N]{};
  static constexpr std::size_t size = N - 1;
  constexpr FixedString(char const (&s)[N]) {
    for (std::size_t i = 0; i < N; ++i)
      data[i] = s[i];
  }
  constexpr char operator[](std::size_t i) const { return data[i]; }
};
template <class X> struct is_Vec : std::false_type {};
template <std::size_t M, class U> struct is_Vec<Vec<M, U>> : std::true_type {};
template <class X>
inline constexpr bool is_Vec_v = is_Vec<std::remove_cvref_t<X>>::value;
template <class X>
inline constexpr bool is_scalar_v =
    std::is_arithmetic_v<std::remove_cvref_t<X>>;
template <class X> struct Vec_size;
template <std::size_t M, class U>
struct Vec_size<Vec<M, U>> : std::integral_constant<std::size_t, M> {};

} // namespace detail

template <std::size_t N, typename T = float>
  requires std::is_arithmetic_v<T>
struct Vec : std::array<T, N> {
private:
  template <class X> static consteval std::size_t extent() {
    if constexpr (detail::is_Vec_v<X>)
      return detail::Vec_size<std::remove_cvref_t<X>>::value;
    else if constexpr (detail::is_scalar_v<X>)
      return 1;
    else
      return 0; // Should be unreachable
  }
  template <class... Args> static consteval std::size_t total_extent() {
    return (extent<Args>() + ... + 0);
  }

public:
  // Constructors
  constexpr Vec() noexcept {
    for (auto &v : *this)
      v = T(0);
  }

  explicit constexpr Vec(T const &s) noexcept {
    for (auto &v : *this)
      v = s;
  }

  template <typename... Args>
    requires((detail::is_scalar_v<Args> || detail::is_Vec_v<Args>) && ...) &&
            (total_extent<Args...>() == N) &&
            (!(sizeof...(Args) == 1 && (detail::is_Vec_v<Args> && ...)))
  constexpr Vec(Args &&...args) noexcept {
    std::size_t i = 0;
    (fill_one(i, std::forward<Args>(args)), ...);
  }

  // Member accesses
  // NOTE: This can (probably) be improved with C++26 reflection in the future.
#define VEC_ACC(component, req, idx)                                           \
  constexpr auto component() noexcept -> T &requires(N >= req) {               \
    return (*this)[idx];                                                       \
  } constexpr auto component() const->T const &                                \
    requires(N >= req)                                                         \
  {                                                                            \
    return (*this)[idx];                                                       \
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

  // RHS operations
  friend constexpr auto operator+(T s, Vec const &v) noexcept -> Vec {
    return v + s;
  }
  friend constexpr auto operator-(T s, Vec const &v) noexcept -> Vec {
    return Vec(s) - v;
  }
  friend constexpr auto operator*(T s, Vec const &v) noexcept -> Vec {
    return v * s;
  }
  friend constexpr auto operator/(T s, Vec const &v) noexcept -> Vec {
    Vec r{};
    for (std::size_t i = 0; i < N; ++i)
      r[i] = s / v[i];
    return r;
  }

  // Members
#define VEC_OP(op)                                                             \
  constexpr auto operator op(Vec const &rhs) const noexcept -> Vec {           \
    Vec result{};                                                              \
    for (std::size_t i = 0; i < N; ++i) {                                      \
      result[i] = (*this)[i] op rhs[i];                                        \
    }                                                                          \
    return result;                                                             \
  }                                                                            \
  constexpr auto operator op(T const &rhs) const noexcept -> Vec {             \
    Vec result{};                                                              \
    for (std::size_t i = 0; i < N; ++i) {                                      \
      result[i] = (*this)[i] op rhs;                                           \
    }                                                                          \
    return result;                                                             \
  }
  VEC_OP(+)
  VEC_OP(-)
  VEC_OP(*)
  VEC_OP(/)
#undef VEC_OP
#define VEC_OP_ASSIGN(sym)                                                     \
  constexpr Vec &operator sym##=(Vec const &rhs) noexcept {                    \
    for (std::size_t i = 0; i < N; ++i)                                        \
      (*this)[i] sym## = rhs[i];                                               \
    return *this;                                                              \
  }                                                                            \
  constexpr Vec &operator sym##=(T const &s) noexcept {                        \
    for (std::size_t i = 0; i < N; ++i)                                        \
      (*this)[i] sym## = s;                                                    \
    return *this;                                                              \
  }
  VEC_OP_ASSIGN(+)
  VEC_OP_ASSIGN(-)
  VEC_OP_ASSIGN(*)
  VEC_OP_ASSIGN(/)
#undef VEC_OP_ASSIGN

  constexpr auto magnitude() const noexcept -> T {
    T total = 0;
    for (auto const &v : *this)
      total += v * v;
    return std::sqrt(total);
  }
  constexpr auto length() const noexcept -> T { return this->magnitude(); }

  constexpr Vec normalized_safe(T eps = eps_default) const noexcept {
    auto m = magnitude();
    return (m > eps) ? (*this) / m : Vec{};
  }
  constexpr Vec normalize_safe(T eps = eps_default) const noexcept {
    return normalized_safe(eps);
  }

  [[nodiscard]] constexpr auto normalized() noexcept -> Vec<N, T> const {
    return (*this) / this->magnitude();
  }
  [[nodiscard]] constexpr auto normalize() noexcept -> Vec<N, T> const {
    return this->normalized();
  }
  [[nodiscard]] constexpr auto unit() noexcept -> Vec<N, T> const {
    return this->normalized();
  }

  [[nodiscard]] constexpr auto dot(Vec<N, T> const &other) const noexcept
      -> T const {
    T res = 0;
    for (std::size_t i = 0; i < N; ++i) {
      res += (*this)[i] * other[i];
    }
    return res;
  }

  static constexpr T eps_default = T(1e-6);
  template <class U = T>
  [[nodiscard]] constexpr auto
  approx_equal(Vec const &rhs, U eps = eps_default) const noexcept {
    using F = std::conditional_t<std::is_floating_point_v<U>, U, double>;
    for (size_t i = 0; i < N; ++i)
      if (std::abs(F((*this)[i] - rhs[i])) > F(eps))
        return false;
    return true;
  }

  template <class U = T> constexpr auto magnitude_promoted() const noexcept {
    using F = std::conditional_t<std::is_floating_point_v<U>, U, double>;
    F s = 0;
    for (auto v : *this)
      s += F(v) * F(v);
    return std::sqrt(s);
  }

  template <typename U = T>
    requires(N == 3)
  constexpr Vec cross(const Vec &r) const noexcept {
    return {(*this)[1] * r[2] - (*this)[2] * r[1],
            (*this)[2] * r[0] - (*this)[0] * r[2],
            (*this)[0] * r[1] - (*this)[1] * r[0]};
  }

  constexpr T distance(Vec const &r) const noexcept {
    return (*this - r).magnitude();
  }

  constexpr Vec project_onto(Vec const &n) const noexcept {
    auto d = this->dot(n);
    auto nn = n.dot(n);
    return (nn ? (d / nn) * n : Vec());
  }

  template <class U>
    requires(std::is_arithmetic_v<U> && N >= 1)
  constexpr explicit(!std::is_convertible_v<U, T>)
      Vec(Vec<N, U> const &other) noexcept {
    for (std::size_t i = 0; i < N; ++i)
      this->operator[](i) = static_cast<T>(other[i]);
  }

  template <class U>
    requires(std::is_arithmetic_v<U> && N >= 1)
  constexpr explicit(!std::is_convertible_v<T, U>)
  operator Vec<N, U>() const noexcept {
    Vec<N, U> r{};
    for (std::size_t i = 0; i < N; ++i)
      r[i] = static_cast<U>((*this)[i]);
    return r;
  }

  template <class U>
    requires(std::is_arithmetic_v<U> && !std::is_same_v<U, T>)
  constexpr Vec &operator=(Vec<N, U> const &rhs) noexcept {
    for (std::size_t i = 0; i < N; ++i)
      (*this)[i] = static_cast<T>(rhs[i]);
    return *this;
  }

private:
  constexpr void fill_one(std::size_t &i, const T &v) noexcept {
    (*this)[i++] = v;
  }
#ifdef SMATH_IMPLICIT_CONVERSIONS
  template <class U>
    requires std::is_arithmetic_v<U> && (!std::is_same_v<U, T>)
  constexpr void fill_one(std::size_t &i, const U &v) noexcept {
    (*this)[i++] = static_cast<T>(v);
  }
  template <std::size_t M, class U>
  constexpr void fill_one(std::size_t &i, const Vec<M, U> &v) noexcept {
    for (std::size_t k = 0; k < M; ++k)
      (*this)[i++] = static_cast<T>(v[k]);
  }
#endif // SMATH_IMPLICIT_CONVERSIONS
  template <std::size_t M>
  constexpr void fill_one(std::size_t &i, const Vec<M, T> &v) noexcept {
    for (std::size_t k = 0; k < M; ++k)
      (*this)[i++] = static_cast<T>(v[k]);
  }
};

template <size_t I, size_t N, class T> constexpr T &get(Vec<N, T> &v) noexcept {
  static_assert(I < N);
  return v[I];
}
template <size_t I, size_t N, class T>
constexpr const T &get(const Vec<N, T> &v) noexcept {
  static_assert(I < N);
  return v[I];
}
template <size_t I, size_t N, class T>
constexpr T &&get(Vec<N, T> &&v) noexcept {
  static_assert(I < N);
  return std::move(v[I]);
}
template <size_t I, size_t N, class T>
constexpr const T &&get(const Vec<N, T> &&v) noexcept {
  static_assert(I < N);
  return std::move(v[I]);
}

template <std::size_t N, typename T = float>
  requires std::is_arithmetic_v<T>
using VecOrScalar = std::conditional_t<N == 1, T, Vec<N, T>>;

namespace detail {

consteval auto char_to_idx(char c) -> std::size_t {
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

constexpr auto is_valid(char c) -> bool {
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
  }
  return false;
}

template <detail::FixedString S, std::size_t N, typename T, std::size_t... I>
constexpr auto swizzle_impl(Vec<N, T> const &v, std::index_sequence<I...>)
    -> VecOrScalar<S.size, T> {
  static_assert(((is_valid(S[I])) && ...), "Invalid swizzle component");
  static_assert(((char_to_idx(S[I]) < N) && ...),
                "Pattern index out of bounds");
  VecOrScalar<S.size, T> out{};
  std::size_t i = 0;
  ((out[i++] = v[char_to_idx(S[I])]), ...);
  return out;
}

template <FixedString S>
concept SwizzleCharsOK = []<std::size_t... I>(std::index_sequence<I...>) {
  return ((is_valid(S[I])) && ...);
}(std::make_index_sequence<S.size>{});

template <FixedString S, std::size_t N>
concept SwizzleInBounds = []<std::size_t... I>(std::index_sequence<I...>) {
  return ((char_to_idx(S[I]) < N) && ...);
}(std::make_index_sequence<S.size>{});

template <FixedString S, std::size_t N>
concept ValidSwizzle =
    (S.size > 0) && SwizzleCharsOK<S> && SwizzleInBounds<S, N>;

} // namespace detail

template <detail::FixedString S, std::size_t N, typename T>
  requires detail::ValidSwizzle<S, N>
constexpr auto swizzle(Vec<N, T> const &v) -> VecOrScalar<S.size, T> {
  return detail::swizzle_impl<S>(v, std::make_index_sequence<S.size>{});
}

using Vec2 = Vec<2>;
using Vec3 = Vec<3>;
using Vec4 = Vec<4>;

using Vec2d = Vec<2, double>;
using Vec3d = Vec<3, double>;
using Vec4d = Vec<4, double>;

} // namespace smath

template <std::size_t N, typename T>
  requires std::formattable<T, char>
struct std::formatter<smath::Vec<N, T>> : std::formatter<T> {
  constexpr auto parse(std::format_parse_context &ctx) {
    return std::formatter<T>::parse(ctx);
  }

  template <typename Ctx>
  auto format(smath::Vec<N, T> const &v, Ctx &ctx) const {
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
template <size_t N, class T>
struct tuple_size<smath::Vec<N, T>> : std::integral_constant<size_t, N> {};

template <size_t I, size_t N, class T>
struct tuple_element<I, smath::Vec<N, T>> {
  static_assert(I < N);
  using type = T;
};
} // namespace std
