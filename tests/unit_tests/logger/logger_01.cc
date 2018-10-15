#include<gtest/gtest.h>
#include <logger.h>

TEST(Logger, default_constructor) {
  Logger l{};
  EXPECT_EQ(l.prefix(), "geotop:");
}

TEST(Logger, push_pop) {
  Logger l{};
  l.push("alberto");
  EXPECT_EQ(l.prefix(), "geotop:alberto:");
  l.pop();
  EXPECT_EQ(l.prefix(), "geotop:");
}

TEST(Logger, geolog) {
  EXPECT_EQ(geolog.prefix(), "geotop:");
}

TEST(ScopedPrefix, basic_test) {
  Logger l{};
  {
    Logger::ScopedPrefix p{"alberto",l};
    EXPECT_EQ(l.prefix(), "geotop:alberto:");
  }
  EXPECT_EQ(l.prefix(), "geotop:");
  {
    Logger::ScopedPrefix p{__func__}; // let us use the function name 
    EXPECT_EQ(geolog.prefix(), "geotop:TestBody:");
  }
  // do the same using the macro
  {
    GEOLOG_PREFIX(__func__); // let us use the function name 
    EXPECT_EQ(geolog.prefix(), "geotop:TestBody:");
  }
  l.pop();
  EXPECT_EQ(l.prefix(), "");
}

TEST(Logger, output_operator){
  testing::internal::CaptureStdout();
  geolog << "ciao" << " alberto" << std::endl
	 << "nuova riga"  << " continua" << std::endl
	 << " queste righe \n sono molto \n piu' complicate" << std::endl
	 <<"con " << std::flush << "jacopo" << std::endl;
#ifndef NDEBUG   
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
#else
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#endif
}

TEST(Logger, file_stream){
  // test that both the channels works
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  geolog.attach_file_stream(std::cerr);
  geolog << "ciao" << " alberto" << std::endl
	 << "nuova riga"  << " continua" << std::endl
	 << " queste righe \n sono molto \n piu' complicate" << std::endl
	 <<"con " << std::flush << "jacopo" << std::endl;
#ifndef NDEBUG   
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
#else
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
#endif
}


TEST(Logger, detach_file_stream){
  // attach and detach a "file" stream
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  geolog.attach_file_stream(std::cerr);
  geolog.detach_file_stream();
  geolog << "ciao" << " alberto" << std::endl
	 << "nuova riga"  << " continua" << std::endl
	 << " queste righe \n sono molto \n piu' complicate" << std::endl
	 <<"con " << std::flush << "jacopo" << std::endl;
#ifndef NDEBUG   
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
#else
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#endif 
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
}

TEST(Logger, iomanip_functions){
  // test the handlig of functions like std::setw, std::setfill
  testing::internal::CaptureStdout();
  std::ostringstream os;
  os << std::setw(10) << 7; 
  geolog << os.str() << std::endl;
  geolog << 1234567890 << std::endl;
#ifndef NDEBUG   
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:         7\ngeotop:1234567890\n");
#else
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#endif   
}

TEST(Logger, for_loop){
  testing::internal::CaptureStdout();
  std::ostringstream os;
  geolog << "entering " << std::endl;
  {
    Logger::ScopedPrefix p {"for"};
    geolog << "vector elements: ";
    for (auto x : {1, 2, 3})
      geolog << x << " ";
    geolog << std::endl;
  }
  geolog << "cycle ended" <<std::endl;
#ifndef NDEBUG   
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:entering \ngeotop:for:vector elements: 1 2 3 \ngeotop:cycle ended\n");
#else
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
#endif   
}




