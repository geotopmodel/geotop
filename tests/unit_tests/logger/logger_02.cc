#include<gtest/gtest.h>
#include <logger.h>
#include <limits>

TEST(Logger, default_depth_levels) {
  Logger l{};
  EXPECT_EQ(l.console_level(), std::numeric_limits<unsigned int>::max());
  EXPECT_EQ(l.file_level(), std::numeric_limits<unsigned int>::max());
  EXPECT_EQ(geolog.console_level(), std::numeric_limits<unsigned int>::max());
  EXPECT_EQ(geolog.file_level(), std::numeric_limits<unsigned int>::max());
  
  EXPECT_EQ(l.prefix_depth_level(), unsigned {1});
  EXPECT_EQ(geolog.prefix_depth_level(), unsigned {1});
}

TEST(Logger, depth_levels) {
  Logger l{};
  l.push("alberto");
  EXPECT_EQ(l.prefix_depth_level(), unsigned {2});

  l.pop();
  EXPECT_EQ(l.prefix_depth_level(), unsigned {1});

  l.pop();
  EXPECT_EQ(l.prefix_depth_level(), unsigned {0});

  EXPECT_ANY_THROW(l.pop());
}


TEST(Logger, ScopedPrefix_depth_levels) {
  Logger l{};
  {
    Logger::ScopedPrefix p{"alberto",l};
    EXPECT_EQ(l.prefix_depth_level(), unsigned {2});
  }
  EXPECT_EQ(l.prefix_depth_level(), unsigned {1});

  {
    Logger::ScopedPrefix p{__func__}; // let us use the function name 
    EXPECT_EQ(geolog.prefix_depth_level(), unsigned {2});
  }
  EXPECT_EQ(geolog.prefix_depth_level(), unsigned {1});
}

TEST(ScopedPrefix, console_depth_level_default) {
  testing::internal::CaptureStdout();
  Logger l{};
  l << "first level"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "second level" << std::endl;
    {
      Logger::ScopedPrefix p{"third",l};
      l << "third level" << std::endl;
    }
    l << "back to second" << std::endl;
  }
  l << "back to first" << std::endl;
  
#ifdef MUTE_GEOLOG
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#else
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:first level\ngeotop:second:second level\ngeotop:second:third:third level\ngeotop:second:back to second\ngeotop:back to first\n");
#endif
  
}


TEST(ScopedPrefix, console_depth_level_0) {
  testing::internal::CaptureStdout();
  Logger l{};

  // no output on screen
  l.set_console_level(0);
  l << "first level"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "second level" << std::endl;
    {
      Logger::ScopedPrefix p{"third",l};
      l << "third level" << std::endl;
    }
    l << "back to second" << std::endl;
  }
  l << "back to first" << std::endl;
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
}

TEST(ScopedPrefix, console_depth_level_1) {
  testing::internal::CaptureStdout();
  Logger l{};

  // only first level
  l.set_console_level(1);
  l << "first level"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "you should not see this" << std::endl;
    {
      Logger::ScopedPrefix p{"third",l};
      l << "nor this" << std::endl;
    }
    l << "nor this one" << std::endl;
  }
  l << "back to first" << std::endl;
  
#ifdef MUTE_GEOLOG
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#else
  EXPECT_EQ("geotop:first level\ngeotop:back to first\n", testing::internal::GetCapturedStdout());
#endif
}

TEST(ScopedPrefix, console_depth_level_2) {
  testing::internal::CaptureStdout();
  Logger l{};

  // first and second level
  l.set_console_level(2);
  l << "first level"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "second level" << std::endl;
    {
      Logger::ScopedPrefix p{"third",l};
      l << "you should not see this" << std::endl;
    }
    l << "back to second" << std::endl;
  }
  l << "back to first" << std::endl;
  
#ifdef MUTE_GEOLOG
    EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#else
    EXPECT_EQ("geotop:first level\ngeotop:second:second level\ngeotop:second:back to second\ngeotop:back to first\n",
	      testing::internal::GetCapturedStdout());
#endif
}


TEST(ScopedPrefix, console_depth_level_changed) {
  testing::internal::CaptureStdout();
  Logger l{};

  l << "first level"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "second level" << std::endl;
    {
      Logger::ScopedPrefix p{"third",l};
      l << "third level" << std::endl;
      l.set_console_level(0);
      // no output from now on
    }
    l << "you should never read this" << std::endl;
  }
  l << "nor this" << std::endl;
  
#ifdef MUTE_GEOLOG
    EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#else
    EXPECT_EQ("geotop:first level\ngeotop:second:second level\ngeotop:second:third:third level\n", testing::internal::GetCapturedStdout());
#endif
}


TEST(ScopedPrefix, console_and_file_levels) {
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  Logger l{};
  l.attach_file_stream(std::cerr);
  l.set_console_level(1);
  l.set_file_level(2);

  l << "this should appear on screen and in file"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "this is only in file" << std::endl;
    {
      Logger::ScopedPrefix p{"third",l};
      l << "you should never see this" << std::endl;
    }
    l << "this is only in file again" << std::endl;
  }
  l << "this should appear on screen and in file again" << std::endl;

#ifdef MUTE_GEOLOG
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
#else
   EXPECT_EQ("geotop:this should appear on screen and in file\ngeotop:this should appear on screen and in file again\n", testing::internal::GetCapturedStdout());
  EXPECT_EQ("geotop:this should appear on screen and in file\ngeotop:second:this is only in file\ngeotop:second:this is only in file again\ngeotop:this should appear on screen and in file again\n",
	    testing::internal::GetCapturedStderr());
#endif
}

