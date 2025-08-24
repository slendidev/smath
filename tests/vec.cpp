#include <format>
#include <string>
#include <type_traits>

#include <gtest/gtest.h>

#include <smath.hpp>

using smath::Vec;
using smath::Vec2;
using smath::Vec3;
using smath::Vec4;

template <class T>
static void ExpectVecNear(const Vec<3, T> &a, const Vec<3, T> &b,
                          T eps = T(1e-6)) {
  for (int i = 0; i < 3; ++i)
    EXPECT_NEAR(double(a[i]), double(b[i]), double(eps));
}

// Constructors and accessors
TEST(Vec, DefaultZero) {
  Vec3 v;
  EXPECT_EQ(v[0], 0.0f);
  EXPECT_EQ(v[1], 0.0f);
  EXPECT_EQ(v[2], 0.0f);
}

TEST(Vec, ScalarFillCtor) {
  Vec4 v{2.0f};
  EXPECT_EQ(v.x(), 2.0f);
  EXPECT_EQ(v.y(), 2.0f);
  EXPECT_EQ(v.z(), 2.0f);
  EXPECT_EQ(v.w(), 2.0f);
}

TEST(Vec, VariadicCtorScalarsAndSubvectors) {
  Vec2 a{1.0f, 2.0f};
  Vec2 b{3.0f, 4.0f};
  Vec4 v{a, b};
  EXPECT_EQ(v.r(), 1.0f);
  EXPECT_EQ(v.g(), 2.0f);
  EXPECT_EQ(v.b(), 3.0f);
  EXPECT_EQ(v.a(), 4.0f);
}

TEST(Vec, NamedAccessorsAliases) {
  Vec3 v{1.0f, 2.0f, 3.0f};
  EXPECT_EQ(v.x(), v.r());
  EXPECT_EQ(v.y(), v.g());
  EXPECT_EQ(v.z(), v.b());
}

// Arithmetic
TEST(Vec, ElementwiseAndScalarOps) {
  Vec3 a{1.0f, 2.0f, 3.0f};
  Vec3 b{4.0f, 5.0f, 6.0f};

  auto s1 = a + b;
  EXPECT_EQ(s1[0], 5.0f);
  EXPECT_EQ(s1[1], 7.0f);
  EXPECT_EQ(s1[2], 9.0f);

  auto s2 = a * 2.0f;
  EXPECT_EQ(s2[0], 2.0f);
  EXPECT_EQ(s2[1], 4.0f);
  EXPECT_EQ(s2[2], 6.0f);

  auto s3 = 2.0f * a; // RHS overloads
  EXPECT_EQ(s3[0], 2.0f);
  EXPECT_EQ(s3[1], 4.0f);
  EXPECT_EQ(s3[2], 6.0f);

  Vec3 c{1.0f, 2.0f, 3.0f};
  c += Vec3{1.0f, 1.0f, 1.0f};
  EXPECT_EQ(c[0], 2.0f);
  EXPECT_EQ(c[1], 3.0f);
  EXPECT_EQ(c[2], 4.0f);

  c *= 2.0f;
  EXPECT_EQ(c[0], 4.0f);
  EXPECT_EQ(c[1], 6.0f);
  EXPECT_EQ(c[2], 8.0f);
}

// Length, dot, cross, normalize
TEST(Vec, MagnitudeAndDot) {
  Vec3 v{3.0f, 4.0f, 12.0f};
  EXPECT_FLOAT_EQ(v.magnitude(), 13.0f);
  EXPECT_FLOAT_EQ(v.length(), 13.0f);

  Vec3 u{1.0f, 0.0f, 2.0f};
  EXPECT_FLOAT_EQ(v.dot(u), 27.0f);
}

TEST(Vec, Cross3D) {
  Vec3 x{1.0f, 0.0f, 0.0f};
  Vec3 y{0.0f, 1.0f, 0.0f};
  auto z = x.cross(y);
  EXPECT_EQ(z[0], 0.0f);
  EXPECT_EQ(z[1], 0.0f);
  EXPECT_EQ(z[2], 1.0f);
}

TEST(Vec, NormalizeAndSafeNormalize) {
  Vec3 v{10.0f, 0.0f, 0.0f};
  auto n = v.normalized();
  auto ns = v.normalized_safe();
  ExpectVecNear(n, Vec3{1.0f, 0.0f, 0.0f});

  Vec3 zero{};
  auto zs = zero.normalized_safe();
  EXPECT_EQ(zs[0], 0.0f);
  EXPECT_EQ(zs[1], 0.0f);
  EXPECT_EQ(zs[2], 0.0f);
}

TEST(Vec, DistanceAndProjection) {
  Vec3 a{1.0f, 2.0f, 3.0f};
  Vec3 b{4.0f, 6.0f, 3.0f};
  EXPECT_FLOAT_EQ(a.distance(b), 5.0f);

  Vec3 n{2.0f, 0.0f, 0.0f};   // onto x-axis scaled
  auto p = a.project_onto(n); // (a·n)/(n·n) * n = (2)/4 * n = 0.5 * n
  ExpectVecNear(p, Vec3{1.0f, 0.0f, 0.0f});
}

// Approx equal
TEST(Vec, ApproxEqual) {
  Vec3 a{1.0f, 2.0f, 3.0f};
  Vec3 b{1.0f + 1e-7f, 2.0f - 1e-7f, 3.0f};
  EXPECT_TRUE(a.approx_equal(b, 1e-6f));
  EXPECT_FALSE(a.approx_equal(b, 1e-9f));
}

// std::get & tuple interop
TEST(Vec, StdGetAndTuple) {
  Vec3 v{7.0f, 8.0f, 9.0f};
  static_assert(std::tuple_size_v<Vec3> == 3);
  static_assert(std::is_same_v<std::tuple_element_t<1, Vec3>, float>);
  EXPECT_EQ(std::get<0>(v), 7.0f);
  EXPECT_EQ(std::get<1>(v), 8.0f);
  EXPECT_EQ(std::get<2>(v), 9.0f);
}

// Swizzle
TEST(Vec, SwizzleBasic) {
  const Vec3 v{1.0f, 2.0f, 3.0f};

  auto yz = smath::swizzle<"yz">(v);
  EXPECT_EQ(yz[0], 2.0f);
  EXPECT_EQ(yz[1], 3.0f);

  auto rxx = smath::swizzle<"xxy">(v);
  EXPECT_EQ(rxx[0], 1.0f);
  EXPECT_EQ(rxx[1], 1.0f);
  EXPECT_EQ(rxx[2], 2.0f);
}

// std::formatter
TEST(Vec, Formatter) {
  smath::Vec<3, int> vi{1, 2, 3};
  std::string s = std::format("{}", vi);
  EXPECT_EQ(s, "{1, 2, 3}");
}

// Conversions
TEST(Vec, ExplicitConversionBetweenScalarTypes) {
  smath::Vec<3, int> vi{1, 2, 3};
  smath::Vec<3, float> vf{vi};
  EXPECT_EQ(vf[0], 1.0f);
  EXPECT_EQ(vf[1], 2.0f);
  EXPECT_EQ(vf[2], 3.0f);

  auto vi2 = static_cast<smath::Vec<3, int>>(vf);
  EXPECT_EQ(vi2[0], 1);
  EXPECT_EQ(vi2[1], 2);
  EXPECT_EQ(vi2[2], 3);
}
