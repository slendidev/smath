# smath

Single-file, header-only linear algebra math library for C++23.

## Features

- Generic `Vec<N, T>` class with useful aliases `Vec2/Vec3/Vec4` and friendly accessors (`x/y/z/w`, `r/g/b/a`). They support approx-equal and tuple/structured bindings.
- `std::format` support.
- Compile-time swizzles via `swizzle<"...">`.
- Generic matrix `Mat` class with useful aliases `Mat2/Mat3/Mat4`.
- `Quaternion<T>` built on `Vec4`.
- Angle helpers `rad/deg/turns` respecting a configurable base unit via the macro `SMATH_ANGLE_UNIT`.
- Optional implicit conversions.
- Packing utilities for normalized RGBA (`pack_unorm4x8`, `unpack_snorm4x8`, etc.).

## License

This library is licensed under the Apache License 2.0. See the (LICENSE.txt)[LICENSE.txt] file for more details.

