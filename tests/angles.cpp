#include <gtest/gtest.h>

#include <smath.hpp>

TEST(AngleReturnRadians, DegInput) {
  EXPECT_NEAR(smath::deg(180.0), std::numbers::pi, 1e-12);
}

TEST(AngleReturnRadians, RadInput) {
  EXPECT_DOUBLE_EQ(smath::rad(std::numbers::pi), std::numbers::pi);
}

TEST(AngleReturnRadians, TurnsInput) {
  EXPECT_NEAR(smath::turns(0.5), std::numbers::pi, 1e-12);
}
