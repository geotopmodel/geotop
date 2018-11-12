INCLUDE(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
  _Pragma(\"GCC diagnostic push\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wextra\\\\\\\"\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wunknown-pragmas\\\\\\\"\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wpragmas\\\\\\\"\")
  int main() { return 0; }
  _Pragma(\"GCC diagnostic pop\")
  "
  COMPILER_HAS_DIAGNOSTIC_PRAGMA)

set(MUTE_GEOLOG OFF CACHE BOOL "Silence the default logger geolog. You may want to turn this to on to gain a bit in terms of performances.")

set(MUTE_GEOTIMER OFF CACHE BOOL "Silence the default timer geotimer. You may want to turn this to on to gain a bit in terms of performances.")


CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/src/geotop/config.cmake.h.in ${CMAKE_BINARY_DIR}/src/geotop/config.h)
# of course we need to inform the compiler where to look for the above header
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/geotop)




set(VERSION ${PROJECT_VERSION})

IF(EXISTS ${CMAKE_SOURCE_DIR}/.git)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
     EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE}  describe --always
     WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
     OUTPUT_VARIABLE "COMMIT"
     ERROR_QUIET
     OUTPUT_STRIP_TRAILING_WHITESPACE)
     MESSAGE( STATUS "Git version: ${GEOtop_BUILD_VERSION}" )
     ELSE()
       SET(COMMIT "Not a git repository")
     ENDIF()
ENDIF()


CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/src/geotop/version.cmake.cc.in ${CMAKE_BINARY_DIR}/src/geotop/version.cc @ONLY)
