
meteoio_prefix = get_option('meteoio_path')

_meteoio_lib = meteoio_prefix+'/lib'

deps += [_cxx.find_library('meteoio', dirs: _meteoio_lib)]
rpaths +=  ':'+_meteoio_lib

include_dirs += [include_directories(meteoio_prefix+'/include')]


