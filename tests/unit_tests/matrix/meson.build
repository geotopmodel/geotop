_tests = ['matrix_01.cc', 'matrix_02.cc']

_src = ['../../../src/geotop/matrix.h']
foreach _t : _tests
    _e = executable(_t, [_t,_src], dependencies: gtest, include_directories: include_dirs)
    test(_t, _e, suite: ['unit_tests','matrix'])
endforeach
