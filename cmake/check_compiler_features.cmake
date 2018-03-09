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
