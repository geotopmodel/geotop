####################
# set additional flag if using intel compiler

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fp_speculation=safe " CACHE STRING "" FORCE)
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-narrowing -Wno-writable-strings -Wno-c++11-narrowing" CACHE STRING "" FORCE)
endif()
