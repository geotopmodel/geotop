main = files('geotop.cc')
_config = files('config.meson.h.in')
src_files += [files('blowingsnow.cc',
                    'channels.cc',
                    'clouds.cc',
                    'deallocate.cc',
                    'energy.balance.cc',
                    'geotop_common.cc',
                    'geotop_common.h',
                    'global_logger.cc',
                    'indices.cc',
                    'input.cc',
                    'inputKeywords.cc',
                    'logger.cc',
                    'meteo.cc',
                    'meteodata.cc',
                    'meteodistr.cc',
                    'output.cc',
                    'output_file.cc',
                    'output_new.cc',
                    'parameters.cc',
                    'PBSM.cc',
                    'pedo.funct.cc',
                    'radiation.cc',
                    'recovering.cc',
                    'snow.cc',
                    'standard_names.cc',
                    'statevar.cc',
                    'tables.cc',
                    'times.cc',
                    'turbulence.cc',
                    'vegetation.cc',
                    'water.balance.cc')]


subdir('Energy')
subdir('Glacier')
subdir('Meteo')
subdir('Snow')
subdir('Soil')
subdir('Water')

include_dirs += [include_directories(['Energy',
                                      'Glacier',
                                      'Meteo',
                                      'Snow',
                                      'Soil',
                                      'Water'])]
