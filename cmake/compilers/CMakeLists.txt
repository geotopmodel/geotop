####################
# set additional flag if using intel compiler

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fp_speculation=safe " CACHE STRING "" FORCE)
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-narrowing -Wno-writable-strings -Wno-c++11-narrowing" CACHE STRING "" FORCE)
endif()


if(WITH_OMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  else()
    error(FATAL_ERROR "Cannot find openMP. Please reconfigure with -DWITH_OMP=OFF")
  endif()
endif()  
