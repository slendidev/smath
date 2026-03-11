#include <type_traits>

#include <gtest/gtest.h>

#include <smath.hpp>

struct ExternalVec3f
{
	float x;
	float y;
	float z;
};

struct ExternalMat2f
{
	float m00;
	float m01;
	float m10;
	float m11;
};

namespace smath
{

template<> struct interop_adapter<ExternalVec3f>
{
	using smath_type = Vec<3, float>;

	static constexpr auto to_smath(ExternalVec3f const &v) -> smath_type
	{
		return { v.x, v.y, v.z };
	}

	static constexpr auto from_smath(smath_type const &v) -> ExternalVec3f
	{
		return { v[0], v[1], v[2] };
	}
};

template<> struct interop_adapter<ExternalMat2f>
{
	using smath_type = Mat<2, 2, float>;

	static constexpr auto to_smath(ExternalMat2f const &m) -> smath_type
	{
		smath_type out {};
		out(0, 0) = m.m00;
		out(0, 1) = m.m01;
		out(1, 0) = m.m10;
		out(1, 1) = m.m11;
		return out;
	}

	static constexpr auto from_smath(smath_type const &m) -> ExternalMat2f
	{
		return { m(0, 0), m(0, 1), m(1, 0), m(1, 1) };
	}
};

} // namespace smath

TEST(Interop, FromExternalVec)
{
	ExternalVec3f ext { 1.0f, 2.0f, 3.0f };
	auto v = smath::from_external(ext);

	static_assert(std::is_same_v<decltype(v), smath::Vec<3, float>>);
	EXPECT_EQ(v[0], 1.0f);
	EXPECT_EQ(v[1], 2.0f);
	EXPECT_EQ(v[2], 3.0f);
}

TEST(Interop, ToExternalVec)
{
	smath::Vec3 v { 4.0f, 5.0f, 6.0f };
	auto ext = smath::to_external<ExternalVec3f>(v);

	EXPECT_EQ(ext.x, 4.0f);
	EXPECT_EQ(ext.y, 5.0f);
	EXPECT_EQ(ext.z, 6.0f);
}

TEST(Interop, RoundtripVec)
{
	ExternalVec3f ext { 7.0f, 8.0f, 9.0f };
	auto v = smath::from_external(ext);
	auto ext2 = smath::to_external<ExternalVec3f>(v);

	EXPECT_EQ(ext2.x, 7.0f);
	EXPECT_EQ(ext2.y, 8.0f);
	EXPECT_EQ(ext2.z, 9.0f);
}

TEST(Interop, MatrixConversion)
{
	ExternalMat2f ext { 1.0f, 2.0f, 3.0f, 4.0f };
	auto m = smath::from_external(ext);

	static_assert(std::is_same_v<decltype(m), smath::Mat<2, 2, float>>);
	EXPECT_EQ(m(0, 0), 1.0f);
	EXPECT_EQ(m(0, 1), 2.0f);
	EXPECT_EQ(m(1, 0), 3.0f);
	EXPECT_EQ(m(1, 1), 4.0f);

	auto ext2 = smath::to_external<ExternalMat2f>(m);
	EXPECT_EQ(ext2.m00, 1.0f);
	EXPECT_EQ(ext2.m01, 2.0f);
	EXPECT_EQ(ext2.m10, 3.0f);
	EXPECT_EQ(ext2.m11, 4.0f);
}
