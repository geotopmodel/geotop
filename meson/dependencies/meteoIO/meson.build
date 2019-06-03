_opt = 'METEOIO'
_dir = get_option(_opt+'_DIR')
_libname = _opt.to_lower()

_code = '''
#include <meteoio/Config.h>
int main()
{}
'''

compiler = meson.get_compiler('cpp')


if _dir == ''

    _compiles = compiler.compiles(_code, name:
                                  _opt+' compile flags check')

    if not _compiles
       error('\n\nCannot find suitable installation for '+_opt+'. Try to specify the prefix with the option -DMETEOIO_DIR')
    endif
    
    _lib = compiler.find_library(_libname,  required: true)

    meteoio = declare_dependency(dependencies:_lib)
else

    _includedir = include_directories(_dir+'/include')
    _compiles = compiler.compiles(_code, name:
                                  _opt+' compile flags check',
                                  include_directories:_includedir)

    if not _compiles
       error('\n\nCannot find suitable installation for '+_opt+' in '+_dir)
    endif

    _libpath = _dir+'/lib'
    _lib = compiler.find_library('meteoio', dirs:_libpath, required: true)

    meteoio = declare_dependency(include_directories: _includedir,
                                 dependencies:_lib)
    
    rpaths = _libpath+':'+rpaths
endif



deps += meteoio


