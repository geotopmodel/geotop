_tests = ['asserts_01.cc']

_src = ['../../../src/geotop/geotop_asserts.h']
foreach _t : _tests
    _e = executable(_t, [_t,_src], dependencies: gtest, include_directories: include_dirs)
    test(_t, _e, suite: ['unit_tests','asserts'])
endforeach
