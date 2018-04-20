#include<gtest/gtest.h>
#include <logger.h>



TEST(Logger, ScopedConsoleLevel) {
  testing::internal::CaptureStdout();
  Logger l{};
  l.set_console_level(0);
  l << "this should not appear"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "this should not appear" << std::endl;
    {
      Logger::ScopedConsoleLevel _l{10,l};
      Logger::ScopedPrefix p{"third",l};
      l << "only this should be printed" << std::endl;
    }
    l << "this should not appear" << std::endl;
  }
  l << "this should not appear" << std::endl;
  EXPECT_EQ("geotop:second:third:only this should be printed\n",testing::internal::GetCapturedStdout());
}

TEST(Logger, ScopedFileLevel) {
  testing::internal::CaptureStderr();
  Logger l{};
  l.attach_file_stream(std::cerr);
  l.set_file_level(0);
  l.set_console_level(0);
  l << "this should not appear"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",l};
    l << "this should not appear" << std::endl;
    {
      Logger::ScopedFileLevel _l{10,l};
      Logger::ScopedPrefix p{"third",l};
      l << "only this should be printed" << std::endl;
    }
    l << "this should not appear" << std::endl;
  }
  l << "this should not appear" << std::endl;
  EXPECT_EQ("geotop:second:third:only this should be printed\n",testing::internal::GetCapturedStderr());
}


// TEST(ScopedPrefix, console_depth_level_0) {
//   testing::internal::CaptureStdout();
//   Logger l{};

//   // no output on screen
//   l.set_console_level(0);
//   l << "first level"<< std::endl;
//   {
//     Logger::ScopedPrefix p{"second",l};
//     l << "second level" << std::endl;
//     {
//       Logger::ScopedPrefix p{"third",l};
//       l << "third level" << std::endl;
//     }
//     l << "back to second" << std::endl;
//   }
//   l << "back to first" << std::endl;
//   EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
// }

// TEST(ScopedPrefix, console_depth_level_1) {
//   testing::internal::CaptureStdout();
//   Logger l{};

//   // only first level
//   l.set_console_level(1);
//   l << "first level"<< std::endl;
//   {
//     Logger::ScopedPrefix p{"second",l};
//     l << "you should not see this" << std::endl;
//     {
//       Logger::ScopedPrefix p{"third",l};
//       l << "nor this" << std::endl;
//     }
//     l << "nor this one" << std::endl;
//   }
//   l << "back to first" << std::endl;
//   EXPECT_EQ("geotop:first level\ngeotop:back to first\n", testing::internal::GetCapturedStdout());
// }

// TEST(ScopedPrefix, console_depth_level_2) {
//   testing::internal::CaptureStdout();
//   Logger l{};

//   // first and second level
//   l.set_console_level(2);
//   l << "first level"<< std::endl;
//   {
//     Logger::ScopedPrefix p{"second",l};
//     l << "second level" << std::endl;
//     {
//       Logger::ScopedPrefix p{"third",l};
//       l << "you should not see this" << std::endl;
//     }
//     l << "back to second" << std::endl;
//   }
//   l << "back to first" << std::endl;
//   EXPECT_EQ("geotop:first level\ngeotop:second:second level\ngeotop:second:back to second\ngeotop:back to first\n",
// 	    testing::internal::GetCapturedStdout());
// }


// TEST(ScopedPrefix, console_depth_level_changed) {
//   testing::internal::CaptureStdout();
//   Logger l{};

//   l << "first level"<< std::endl;
//   {
//     Logger::ScopedPrefix p{"second",l};
//     l << "second level" << std::endl;
//     {
//       Logger::ScopedPrefix p{"third",l};
//       l << "third level" << std::endl;
//       l.set_console_level(0);
//       // no output from now on
//     }
//     l << "you should never read this" << std::endl;
//   }
//   l << "nor this" << std::endl;
//   EXPECT_EQ("geotop:first level\ngeotop:second:second level\ngeotop:second:third:third level\n", testing::internal::GetCapturedStdout());
// }


// TEST(ScopedPrefix, console_and_file_levels) {
//   testing::internal::CaptureStdout();
//   testing::internal::CaptureStderr();
//   Logger l{};
//   l.attach_file_stream(std::cerr);
//   l.set_console_level(1);
//   l.set_file_level(2);

//   l << "this should appear on screen and in file"<< std::endl;
//   {
//     Logger::ScopedPrefix p{"second",l};
//     l << "this is only in file" << std::endl;
//     {
//       Logger::ScopedPrefix p{"third",l};
//       l << "you should never see this" << std::endl;
//     }
//     l << "this is only in file again" << std::endl;
//   }
//   l << "this should appear on screen and in file again" << std::endl;

//   EXPECT_EQ("geotop:this should appear on screen and in file\ngeotop:this should appear on screen and in file again\n", testing::internal::GetCapturedStdout());
//   EXPECT_EQ("geotop:this should appear on screen and in file\ngeotop:second:this is only in file\ngeotop:second:this is only in file again\ngeotop:this should appear on screen and in file again\n", testing::internal::GetCapturedStderr());
// }

