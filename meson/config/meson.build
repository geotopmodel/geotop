include_dirs = [include_dirs, include_directories('.')]

code = '''_Pragma("GCC diagnostic push")
  _Pragma("GCC diagnostic ignored \"-Wextra\"")
  _Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")
  _Pragma("GCC diagnostic ignored \"-Wpragmas\"")
  int main() { return 0; }
  _Pragma("GCC diagnostic pop")
'''

_pragma = _cxx.compiles(code, name: 'pragma check')

_cdata = configuration_data()

configure_file(output: 'config.h',
               input: _config,
               configuration: _cdata)

_cdata = configuration_data()



_project_git_version = '0'
# test if we have git and if it we have it, we query the git version
# of our project
_git = run_command('which','git')
if _git.returncode() == 0
    r = run_command('git', 'describe', '--tags', '--always', '--abbrev=8')
    _project_git_version = r.stdout().strip()
    message('Git version: '+_project_git_version)
endif

_cdata.set_quoted('GEOtop_BUILD_VERSION', _project_git_version)
_cdata.set_quoted('PACKAGE_STRING', meson.project_name()+' '+meson.project_version())

configure_file(output: 'version.h',
               configuration: _cdata)
