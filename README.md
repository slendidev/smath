<p align="center">
	<img alt="smath logo" width="100%" src="https://paste.slendi.dev/LcCPsJ.svg" />
	<br><br>
	<a href="https://github.com/slendidev/smath/actions/workflows/build.yml">
		<img src="https://github.com/slendidev/smath/actions/workflows/build.yml/badge.svg" alt="Build">
	</a>
	<a href="LICENSE.txt">
		<img src="https://img.shields.io/badge/license-Apache%202.0-blue.svg" alt="License: Apache-2.0">
	</a>
	<a href="https://en.cppreference.com/w/cpp/23">
		<img src="https://img.shields.io/badge/C%2B%2B-23-00599C.svg" alt="C++23">
	</a>
</p>

## Features

- Generic `Vec<N, T>` class with useful aliases `Vec2/Vec3/Vec4` and friendly accessors (`x/y/z/w`, `r/g/b/a`). They support approx-equal and tuple/structured bindings.
- `std::format` support.
- Compile-time swizzles via `swizzle<"...">`.
- Generic matrix `Mat` class with useful aliases `Mat2/Mat3/Mat4`.
- `Quaternion<T>` built on `Vec4`.
- Angle helpers `rad/deg/turns` respecting a configurable base unit via the macro `SMATH_ANGLE_UNIT`.
- Optional implicit conversions.
- Packing utilities for normalized RGBA (`pack_unorm4x8`, `unpack_snorm4x8`, etc.).
- C++20 modules support.
- Built-in interop adapter and optional interop headers for various libraries.

## Quick Start

Create `main.cpp`:

```cpp
#include <print>
#include <smath/smath.hpp>

int main() {
  using namespace smath;

  Vec3 a{1.0f, 2.0f, 3.0f};
  Vec3 b{3.0f, 2.0f, 1.0f};

  std::println("a + b = {}", a + b);
  std::println("dot(a, b) = {}", a.dot(b));
  std::println("normalized(a) = {}", a.normalized());
}
```

Build and run:

```bash
g++ -std=c++23 -Iinclude main.cpp -o quickstart
./quickstart
```

## Interop Headers

smath ships optional interop adapters under `include/smath/interop/`.

When installing with CMake, interop headers are gated behind:

- `-DSMATH_INSTALL_INTEROP_HEADERS=ON`

All of them plug into the same generic API defined in `smath.hpp`:

```cpp
auto v = smath::from_external(external_value);
auto ext = smath::to_external<ExternalType>(smath_value);
```

## License

This library is licensed under the Apache License 2.0. See the [LICENSE.txt](LICENSE.txt) file for more details.
