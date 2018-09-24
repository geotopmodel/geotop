_tests = ['tensor-debug_01.cc']

_src = ['../../../src/geotop/tensor.h',
	'../../../src/libraries/fluidturtle/turtle.h',
	'../../../src/libraries/fluidturtle/tensors3D.cc',
	'../../../src/libraries/fluidturtle/error.cc']

foreach _t : _tests
    _e = executable(_t, [_t,_src], dependencies: gtest, include_directories: include_dirs)
    test(_t, _e, suite: ['unit_tests','tensor-debug'])
endforeach
