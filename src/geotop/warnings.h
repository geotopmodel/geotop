#ifdef COMPILER_HAS_DIAGNOSTIC_PRAGMA
#define DISABLE_WARNINGS                                                \
  _Pragma("GCC diagnostic push")                                          \
  _Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")                 \
  _Pragma("GCC diagnostic ignored \"-Wpragmas\"")                         \
  _Pragma("GCC diagnostic ignored \"-Wunknown-warning-option\"")          \
  _Pragma("GCC diagnostic ignored \"-Wunknown-warning\"")                 \
  _Pragma("GCC diagnostic ignored \"-Wextra\"")                           \
  _Pragma("GCC diagnostic ignored \"-Waddress-of-packed-member\"")        \
  _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")         \
  _Pragma("GCC diagnostic ignored \"-Wexpansion-to-defined\"")            \
  _Pragma("GCC diagnostic ignored \"-Wexpansion-to-defined\"")            \
  _Pragma("GCC diagnostic ignored \"-Wignored-attributes\"")              \
  _Pragma("GCC diagnostic ignored \"-Wignored-qualifiers\"")              \
  _Pragma("GCC diagnostic ignored \"-Wimplicit-fallthrough\"")            \
  _Pragma("GCC diagnostic ignored \"-Winfinite-recursion\"")              \
  _Pragma("GCC diagnostic ignored \"-Wmisleading-indentation\"")          \
  _Pragma("GCC diagnostic ignored \"-Wmissing-field-initializers\"")      \
  _Pragma("GCC diagnostic ignored \"-Wnon-virtual-dtor\"")                \
  _Pragma("GCC diagnostic ignored \"-Woverflow\"")                        \
  _Pragma("GCC diagnostic ignored \"-Woverloaded-virtual\"")              \
  _Pragma("GCC diagnostic ignored \"-Wpedantic\"")                        \
  _Pragma("GCC diagnostic ignored \"-Wstrict-aliasing\"")                 \
  _Pragma("GCC diagnostic ignored \"-Wtautological-constant-out-of-range-compare\"") \
  _Pragma("GCC diagnostic ignored \"-Wtype-limits\"")                     \
  _Pragma("GCC diagnostic ignored \"-Wundef\"")                           \
  _Pragma("GCC diagnostic ignored \"-Wunused-but-set-parameter\"")        \
  _Pragma("GCC diagnostic ignored \"-Wunused-but-set-variable\"")         \
  _Pragma("GCC diagnostic ignored \"-Wunused-function\"")                 \
  _Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")                \
  _Pragma("GCC diagnostic ignored \"-Wunused-private-field\"")            \
  _Pragma("GCC diagnostic ignored \"-Wunused-variable\"")                 \
  _Pragma("GCC diagnostic ignored \"-Wpragmas\"")

#define ENABLE_WARNINGS \
  _Pragma("GCC diagnostic pop")

#else

#define DISABLE_WARNINGS
#define ENABLE_WARNINGS

#endif
