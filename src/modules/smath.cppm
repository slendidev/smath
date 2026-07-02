module;

#include "smath/smath.hpp"

export module smath;

export namespace smath
{
using ::smath::Mat;
using ::smath::Quaternion;
using ::smath::Vec;
using ::smath::VecOrScalar;

using ::smath::Vec2;
using ::smath::Vec3;
using ::smath::Vec4;

using ::smath::Vec2d;
using ::smath::Vec3d;
using ::smath::Vec4d;

using ::smath::Mat2;
using ::smath::Mat3;
using ::smath::Mat4;

using ::smath::Mat2d;
using ::smath::Mat3d;
using ::smath::Mat4d;

using ::smath::swizzle;

using ::smath::deg;
using ::smath::rad;
using ::smath::turns;

using ::smath::pack_snorm4x8;
using ::smath::pack_unorm4x8;
using ::smath::unpack_snorm4x8;
using ::smath::unpack_unorm4x8;

using ::smath::operator*;

using ::smath::rotate;
using ::smath::scale;
using ::smath::translate;

using ::smath::shear_x;
using ::smath::shear_y;
using ::smath::shear_z;

using ::smath::matrix_infinite_perspective;
using ::smath::matrix_look_at;
using ::smath::matrix_ortho3d;
using ::smath::matrix_perspective;
} // namespace smath
