# The meson build system for the geotop project
# Author: Alberto Sartori (alberto.sartori86@gmail.com) 2018

project('geotop','cpp',
        default_options : ['c++std=c11',
			   'buildtype=release',
			   'warning_level=3' ],
        version: '2.0.0')


include_dirs = include_directories ([ 'src/geotop',
                                      'src/libraries/ascii',
                                      'src/libraries/fluidturtle',
                                      'src/libraries/geomorphology',
                                      'src/libraries/math'])
src_files = ([
        'src/geotop/blowingsnow.cc',
        'src/geotop/channels.cc',
        'src/geotop/clouds.cc',
        'src/geotop/deallocate.cc',
        'src/geotop/energy.balance.cc',
        'src/geotop/indices.cc',
        'src/geotop/input.cc',
        'src/geotop/meteo.cc',
        'src/geotop/meteodata.cc',
        'src/geotop/meteodistr.cc',
        'src/geotop/output.cc',
        'src/geotop/parameters.cc',
        'src/geotop/PBSM.cc',
        'src/geotop/pedo.funct.cc',
        'src/geotop/radiation.cc',
        'src/geotop/recovering.cc',
        'src/geotop/snow.cc',
        'src/geotop/tables.cc',
        'src/geotop/times.cc',
        'src/geotop/turbulence.cc',
        'src/geotop/vegetation.cc',
        'src/geotop/water.balance.cc',
        'src/libraries/ascii/import_ascii.cc',
        'src/libraries/ascii/init.cc',
        'src/libraries/ascii/rw_maps.cc',
        'src/libraries/ascii/tabs.cc',
        'src/libraries/ascii/write_ascii.cc',
        'src/libraries/fluidturtle/alloc.cc',
        'src/libraries/fluidturtle/datamanipulation.cc',
        'src/libraries/fluidturtle/error.cc',
        'src/libraries/fluidturtle/linearalgebra.cc',
        'src/libraries/fluidturtle/list.cc',
        'src/libraries/fluidturtle/probability.cc',
        'src/libraries/fluidturtle/random.cc',
        'src/libraries/fluidturtle/statistics.cc',
        'src/libraries/fluidturtle/t_io.cc',
        'src/libraries/fluidturtle/tensors3D.cc',
        'src/libraries/fluidturtle/utilities.cc',
        'src/libraries/fluidturtle/write_dem.cc',
        'src/libraries/geomorphology/dtm_resolution.cc',
        'src/libraries/geomorphology/geomorphology.0875.cc',
        'src/libraries/geomorphology/geomorphology.cc',
        'src/libraries/geomorphology/networks.cc',
        'src/libraries/geomorphology/shadows.cc',
        'src/libraries/math/sparse_matrix.cc',
        'src/libraries/math/util_math.cc'])


main = files('src/geotop/geotop.cc')
# m = meson.get_compiler('c').find_library('m')

geotop = executable('geotop-2.0.0', main, src_files,
                    include_directories: include_dirs)
                    # dependencies:m)

