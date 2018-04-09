_prefix = meson.current_source_dir()
# list all the tests with the following format
#     folder/test_name
# each test will have the following properties:
# name : 'folder/test_name'
# suite : ['folder', 'test_name', name]

# you can run a specific test with its name attribute
#     meson test name
# i.e.
#     meson test folder/test_name

# or a set of tests by means of a suite 
#     meson test --suite geotop:folder
#     meson test --suite geotop:test_name

_tests = ['1D/B2_BeG_017',
          '1D/Bro',
          '1D/Calabria',
          '1D/ColdelaPorte',
          '1D/CostantMeteo',
          '1D/Jungfraujoch',
          '1D/Matsch_B2_Ref_007',
          '1D/Matsch_P2_Ref_007',
          '1D/PureDrainage',
          '1D/PureDrainageFaked',
          '1D/PureDrainageRainy',
          '1D/PureDrainageRainySlope',
          '3D/Borden05m',
          '3D/example',
          '3D/hillslope01',
          '3D/hillslope02_superslab',
          '3D/Mazia',
          '3D/no_reflection',
          '3D/onepoint_hydrostatic',
          '3D/panola',
          '3D/panola_25pixel',
          '3D/panola_25pixel_nobed',
          '3D/panola_25pixel_nobed_hydrostatic',
          '3D/prealpiC',
          '3D/PSQL_test',
          '3D/rendena',
          '3D/small_example',
          '3D/small_example-channel',
          '3D/small_example-onlyEnergy',
          '3D/snow_dstr_SENSITIVITY',
          '3D/Vshape',
          '3D/WG1_2.0_001']

foreach _t : _tests
    _test_prefix = _prefix+'/'+_t
   test(_t, geotop, args:[_test_prefix], suite: _t.split('/'),
        timeout: 900)
endforeach

# setup the test_runner tests. They have an additional suite entry,
# namely 'test_runner'

_conf_runner = configuration_data()
_conf_runner.set('GEOTOP_PROGRAM_PATH', meson.build_root()+'/geotop')

_meteoio_status = 'off'
if get_option('meteoio_output')
    _meteoio_status = 'on'   
endif
_conf_runner.set('METEOIOPREFIX', _meteoio_status)

_conf_runner.set('GEOTOP_TEST_DIR', _prefix)

configure_file(output: 'test_runner.py',
               input: 'test_runner.py.in',
               configuration: _conf_runner)

py2 = find_program('python2')
tr = meson.current_build_dir()+'/test_runner.py'

foreach _t : _tests
    _test_prefix = _prefix+'/'+_t
    test(_t+'.test_runner', py2, args:[tr],
         workdir: _test_prefix, is_parallel: false, suite: _t.split('/'),
         timeout: 900)
endforeach
