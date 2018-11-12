#include <gtest/gtest.h>
#include <timer.h>

TEST(ClockMeasurements, default_constructor){
  ClockMeasurements c;
  EXPECT_EQ(0u, c.number_of_calls);
  EXPECT_EQ(0u, c.elapsed_time);
}
