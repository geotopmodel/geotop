#include<gtest/gtest.h>
#include <logger.h>


TEST(Logger, default_constructor) {
  Logger l{};
  EXPECT_EQ(l.pop(), "geotop:");
}

TEST(Logger, push_pop) {
  Logger l{};
  l.push("alberto");
  EXPECT_EQ(l.pop(), "geotop:alberto:");
  EXPECT_EQ(l.pop(), "geotop:");
}

TEST(Logger, geolog) {
  EXPECT_EQ(geolog.prefix(), "geotop:");
}

TEST(Prefix, basic_test) {
  Logger l{};
  {
    Logger::Prefix p{"alberto",l};
    EXPECT_EQ(l.prefix(), "geotop:alberto:");
  }
  EXPECT_EQ(l.prefix(), "geotop:");
}

TEST(Logger, output_operator){
  testing::internal::CaptureStdout();
  geolog << "ciao" << " alberto" << std::endl
	 << "nuova riga"  << " continua" << std::endl
	 << " queste righe \n sono molto \n piu' complicate" << std::endl
	 <<"con " << std::flush << "jacopo" << std::endl;
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
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
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
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
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "geotop:ciao alberto\ngeotop:nuova riga continua\ngeotop: queste righe \n sono molto \n piu' complicate\ngeotop:con jacopo\n");
  EXPECT_EQ(testing::internal::GetCapturedStderr(), "");
}


