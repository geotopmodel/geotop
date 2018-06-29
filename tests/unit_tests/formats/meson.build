_tests = ['date.cc', 'parameter.cc']

foreach _t : _tests
    _e = executable(_t, [_t], dependencies: gtest, include_directories: include_dirs)
    test(_t, _e, suite: ['unit_tests','formats'])
endforeach
