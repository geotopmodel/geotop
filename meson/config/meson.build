include_dirs += include_directories('.')

# build features
_cdata = configuration_data()

# user defined options
foreach _opt : ['MUTE_GEOLOG', 'MUTE_GEOTIMER']
    _cdata.set(_opt, get_option(_opt))
endforeach


# check if compiler supports pragmas
code = '''_Pragma("GCC diagnostic push")
  _Pragma("GCC diagnostic ignored \"-Wextra\"")
  _Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")
  _Pragma("GCC diagnostic ignored \"-Wpragmas\"")
  int main() { return 0; }
  _Pragma("GCC diagnostic pop")
'''

_pragma = _cxx.compiles(code, name: 'pragma check')

_cdata.set('COMPILER_HAS_DIAGNOSTIC_PRAGMA', _pragma)

configure_file(output: 'config.h',
               input: _config,
               configuration: _cdata)


# take care of reporting proper version number and commit
_cdata = configuration_data()

_cdata.set('VERSION',meson.project_version() )

_v = configure_file(input: _vcs,
               output:'version.cc.1',
               configuration: _cdata)

tag = vcs_tag( command: ['git', 'describe', '--always'],
               input: _v,
               output: 'version.cc',
               replace_string: 'XXX',
               fallback:'Not a git repository')
