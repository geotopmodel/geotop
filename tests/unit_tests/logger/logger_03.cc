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
  
#ifndef NDEBUG
  EXPECT_EQ("geotop:second:third:only this should be printed\n",testing::internal::GetCapturedStdout());
#else
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
#endif
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
  
#ifndef NDEBUG  
  EXPECT_EQ("geotop:second:third:only this should be printed\n",testing::internal::GetCapturedStderr());
#else
  EXPECT_EQ("", testing::internal::GetCapturedStderr());
#endif
}

TEST(Logger, ScopedLevels){
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  Logger log{};
  log.attach_file_stream(std::cerr);
  //everything should be suppressed
  Logger::ScopedLevels _l{0,log};
  log << "this should not appear"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",log};
    log << "this should not appear" << std::endl;
    {
      Logger::ScopedPrefix p{"third",log};
      log << "this should not appear" << std::endl;
    }
    log << "this should not appear" << std::endl;
  }
  log << "this should not appear" << std::endl;

  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  EXPECT_EQ("", testing::internal::GetCapturedStderr());
}

TEST(Logger, ScopedLevels_mixed_01){
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  Logger log{};

  Logger::ScopedLevels _l{2,0,log};
  log << "this should appear on stdout"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",log};
    log << "this should appear on stdout" << std::endl;
    {
      Logger::ScopedPrefix p{"third",log};
      log << "this should not appear" << std::endl;
    }
    log << "this should appear on stdout" << std::endl;
  }
  log << "this should appear on stdout" << std::endl;

#ifndef NDEBUG  
  EXPECT_EQ("geotop:this should appear on stdout\ngeotop:second:this should appear on stdout\ngeotop:second:this should appear on stdout\ngeotop:this should appear on stdout\n", testing::internal::GetCapturedStdout());
#else
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
#endif
  EXPECT_EQ("", testing::internal::GetCapturedStderr());
  
}

TEST(Logger, ScopedLevels_mixed_02){
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  Logger log{};

  log.attach_file_stream(std::cerr);
  
  Logger::ScopedLevels _l{2,1,log};
  log << "this should appear on stdout and stderr"<< std::endl;
  {
    Logger::ScopedPrefix p{"second",log};
    log << "this should appear on stdout" << std::endl;
    {
      Logger::ScopedPrefix p{"third",log};
      log << "this should not appear" << std::endl;
    }
    log << "this should appear on stdout" << std::endl;
  }
  log << "this should appear on stdout and stderr" << std::endl;

#ifndef NDEBUG  
  EXPECT_EQ("geotop:this should appear on stdout and stderr\ngeotop:second:this should appear on stdout\ngeotop:second:this should appear on stdout\ngeotop:this should appear on stdout and stderr\n", testing::internal::GetCapturedStdout());
  EXPECT_EQ("geotop:this should appear on stdout and stderr\ngeotop:this should appear on stdout and stderr\n", testing::internal::GetCapturedStderr());
#else
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  EXPECT_EQ("", testing::internal::GetCapturedStderr());
#endif
}


TEST(ScopedLevels, ScopedLevels_mixed_03){
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  Logger log{};
  log.attach_file_stream(std::cerr);
  //everything should be suppressed
  Logger::ScopedLevels _l{0,log};
  log << "this should not appear"<< std::endl;
  {
    Logger::ScopedLevels _sl{10,0,log};
    Logger::ScopedPrefix p{"second",log};
    log << "this should appear on stdout" << std::endl;
    {
      Logger::ScopedLevels _ssl{0,10,log};
      Logger::ScopedPrefix p{"third",log};
      log << "this should appear on stderr" << std::endl;
    }
    log << "this should appear on stdout again" << std::endl;
  }
  log << "this should not appear" << std::endl;

#ifndef NDEBUG  
  EXPECT_EQ("geotop:second:this should appear on stdout\ngeotop:second:this should appear on stdout again\n", testing::internal::GetCapturedStdout());
  EXPECT_EQ("geotop:second:third:this should appear on stderr\n", testing::internal::GetCapturedStderr());
#else
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  EXPECT_EQ("", testing::internal::GetCapturedStderr());
#endif
}


TEST(ScopedLevels, mixed_geolog){
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  geolog.attach_file_stream(std::cerr);
  //everything should be suppressed
  Logger::ScopedLevels _l{0};
  geolog << "this should not appear"<< std::endl;
  {
    Logger::ScopedLevels _sl{10,0};
    Logger::ScopedPrefix p{"second"};
    geolog << "this should appear on stdout" << std::endl;
    {
      Logger::ScopedLevels _ssl{0,10};
      Logger::ScopedPrefix p{"third"};
      geolog << "this should appear on stderr" << std::endl;
    }
    geolog << "this should appear on stdout again" << std::endl;
  }
  geolog << "this should not appear" << std::endl;
  
#ifndef NDEBUG  
  EXPECT_EQ("geotop:second:this should appear on stdout\ngeotop:second:this should appear on stdout again\n", testing::internal::GetCapturedStdout());
  EXPECT_EQ("geotop:second:third:this should appear on stderr\n", testing::internal::GetCapturedStderr());
#else
  EXPECT_EQ("", testing::internal::GetCapturedStdout());
  EXPECT_EQ("", testing::internal::GetCapturedStderr());
#endif
}
