#include<gtest/gtest.h>
#include <iostream>
#include <iomanip>

TEST(Date, format){
  // we print the dates in the format
  // dd/mm/yyyy hh:mm
  // here we test that the output is as expected
  // and we properly use the functions of iomanip
  testing::internal::CaptureStdout();
  int i=6;
  std::cout << std::setfill('0') << std::setw(2) << i;
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "06");
}

