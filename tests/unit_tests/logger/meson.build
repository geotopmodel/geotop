_tests = ['logger_01.cc',
          'logger_02.cc']
          # 'logger_03.cc']

_src = ['../../../src/geotop/logger.cc']

foreach _t : _tests
    _e = executable(_t, _t,_src, dependencies: gtest, include_directories: include_dirs)
    test(_t, _e, suite: ['unit_tests','logger'])
endforeach
