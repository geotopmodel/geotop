PROJECT(utest)

IF(COMMAND CMAKE_POLICY)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

ADD_EXECUTABLE(utest
	maintest.cpp
	test_input_keywords.cpp
	)

TARGET_LINK_LIBRARIES(utest
	
	${INPUT_KEYWORDS_LIBRARIES}
	)

#LIBS += -lpthread
#LIBS += -lboost_unit_test_framework-mt
#LIBS += -lboost_thread-mt
#LIBS += -lboost_system-mt
