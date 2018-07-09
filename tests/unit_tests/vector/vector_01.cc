#include <gtest/gtest.h>
#include <vector.h>


TEST(Vector, constructors){
  Vector<int> v{3};
  EXPECT_EQ(v.size(), std::size_t{3});
  
  Vector<int> v0{3,0};
  EXPECT_EQ(std::size_t{4}, v0.size());

  Vector<int> vcopied{v0};
  EXPECT_EQ(vcopied.size(), std::size_t{4});
}

TEST(Vector, range_for){
  Vector<double> v{3};
  double c{1.0};
  v[0] = -9999;
  v[3] = -9999;

  for (auto &x : v)
    x = ++c;

  EXPECT_DOUBLE_EQ(v[0], -9999);
  EXPECT_DOUBLE_EQ(v[1], 2);
  EXPECT_DOUBLE_EQ(v[2], 3);
  EXPECT_DOUBLE_EQ(v[3], 4);

}

TEST(Vector, range_for_zero){
  Vector<double> v{3,0};
  double c{0.0};
  v[0] = -9999;
  v[3] = -9999;
  
  for (auto &x : v)
    x = ++c;

  EXPECT_DOUBLE_EQ(v[0], 1);
  EXPECT_DOUBLE_EQ(v[1], 2);
  EXPECT_DOUBLE_EQ(v[2], 3);
  EXPECT_DOUBLE_EQ(v[3], 4);

}

TEST(Vector, set_value){
  Vector<double> v{3,0};
  double c{0.0};
  
  for (auto &x : v)
    x = ++c;

  v = -9999.;

  EXPECT_DOUBLE_EQ(v[0], -9999.);
  EXPECT_DOUBLE_EQ(v[1], -9999.);
  EXPECT_DOUBLE_EQ(v[2], -9999.);
  EXPECT_DOUBLE_EQ(v[3], -9999.);
}



TEST(Vector, copy_semantic){
  Vector<double> v{3,0};
  double c{0.0};
  v[0] = -9999.;
  v[3] = -9999.;
  
  for (auto &x : v)
    x = ++c;
  
  Vector<double> v1{v};

  v = -9999.0; // change v to check the vectors are not linked
  EXPECT_DOUBLE_EQ(v1[0], 1);
  EXPECT_DOUBLE_EQ(v1[1], 2);
  EXPECT_DOUBLE_EQ(v1[2], 3);
  EXPECT_DOUBLE_EQ(v1[3], 4);

  Vector<double> v2{1};
  v2[1] = -9999;

  for (auto &x : v)
    x = ++c;

  v2=v1;
  
  v = -9999.0; // change v1 to check the vectors are not linked

  ASSERT_EQ(v1.size(), v2.size());
  EXPECT_DOUBLE_EQ(v2[0], 1);
  EXPECT_DOUBLE_EQ(v2[1], 2);
  EXPECT_DOUBLE_EQ(v2[2], 3);
  EXPECT_DOUBLE_EQ(v2[3], 4);
  
}


TEST(Vector, out_of_range){
  Vector<double> v{3};
  double c{1.0};
  
  EXPECT_NO_THROW(v.at(1));
  EXPECT_NO_THROW(v.at(2));
  EXPECT_NO_THROW(v.at(3));
  
  EXPECT_ANY_THROW(v.at(0));
  EXPECT_ANY_THROW(v.at(4));
  EXPECT_ANY_THROW(v.at(-3));
  EXPECT_ANY_THROW(v.at(5000));


  EXPECT_NO_THROW(v(1));
  EXPECT_NO_THROW(v(2));
  EXPECT_NO_THROW(v(3));

#ifndef NDEBUG
  EXPECT_ANY_THROW(v(0));
  EXPECT_ANY_THROW(v(4));
  EXPECT_ANY_THROW(v(-3));
  EXPECT_ANY_THROW(v(5000));
#else
  EXPECT_NO_THROW(v(0));
  EXPECT_NO_THROW(v(4));
  EXPECT_NO_THROW(v(-3));
  EXPECT_NO_THROW(v(5000));
#endif
  
}

TEST(Vector, out_of_range_zero){
  Vector<double> v{3,0};
  double c{1.0};

  EXPECT_NO_THROW(v.at(0));
  EXPECT_NO_THROW(v.at(1));
  EXPECT_NO_THROW(v.at(2));
  EXPECT_NO_THROW(v.at(3));
  
  EXPECT_ANY_THROW(v.at(-3));
  EXPECT_ANY_THROW(v.at(4));
  EXPECT_ANY_THROW(v.at(5000));
  
}

TEST(Vector, summation){
  Vector<double> v_0{3};
  
  Vector<double> v_1{3,0};

  // vector lenght mismatch
#ifndef NDEBUG
  EXPECT_ANY_THROW(v_1 += v_0);
#else
  EXPECT_NO_THROW(v_1 += v_0);
#endif
  
  Vector<double> v_2{3,0};
  double c{0.0};  
  for (auto &x : v_1)
    x = ++c;

  for (auto &x : v_2)
    x = ++c;

  EXPECT_DOUBLE_EQ(v_2[0], 5.);
  EXPECT_DOUBLE_EQ(v_2[1], 6.);
  EXPECT_DOUBLE_EQ(v_2[2], 7.);
  EXPECT_DOUBLE_EQ(v_2[3], 8.);

  EXPECT_NO_THROW(v_2 += v_1);
  
  EXPECT_DOUBLE_EQ(v_1[0], 1.);
  EXPECT_DOUBLE_EQ(v_1[1], 2.);
  EXPECT_DOUBLE_EQ(v_1[2], 3.);
  EXPECT_DOUBLE_EQ(v_1[3], 4.);

  EXPECT_DOUBLE_EQ(v_2[0], 6.);
  EXPECT_DOUBLE_EQ(v_2[1], 8.);
  EXPECT_DOUBLE_EQ(v_2[2], 10.);
  EXPECT_DOUBLE_EQ(v_2[3], 12.);
}
