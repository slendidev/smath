module;

#include "smath.hpp"

export module smath;

export namespace smath {
export using ::smath::Vec;
export using ::smath::VecOrScalar;
export using ::smath::Quaternion;
export using ::smath::Mat;

export using ::smath::Vec2;
export using ::smath::Vec3;
export using ::smath::Vec4;

export using ::smath::Vec2d;
export using ::smath::Vec3d;
export using ::smath::Vec4d;

export using ::smath::Mat2;
export using ::smath::Mat3;
export using ::smath::Mat4;

export using ::smath::Mat2d;
export using ::smath::Mat3d;
export using ::smath::Mat4d;

export using ::smath::swizzle;

export using ::smath::deg;
export using ::smath::rad;
export using ::smath::turns;

export using ::smath::pack_unorm4x8;
export using ::smath::pack_snorm4x8;
export using ::smath::unpack_unorm4x8;
export using ::smath::unpack_snorm4x8;

export using ::smath::operator*;

export using ::smath::translate;
export using ::smath::rotate;
export using ::smath::scale;

export using ::smath::shear_x;
export using ::smath::shear_y;
export using ::smath::shear_z;

export using ::smath::matrix_ortho3d;
export using ::smath::matrix_perspective;
export using ::smath::matrix_infinite_perspective;
export using ::smath::matrix_look_at;
}
