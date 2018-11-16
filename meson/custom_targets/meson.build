# define indent target for indenting all source files of geotop

_source_dir = meson.source_root()
_build_dir  = meson.build_root()

_indent = find_program([_source_dir+'/utils/meson_indent'], required: false )

run_target('indent',
           command: _indent)

_run_clang_tidy = find_program(_source_dir+'/utils/run_clang_tidy', required: false)

if _run_clang_tidy.found()
    run_target('tidy',command:
               [_run_clang_tidy,_build_dir, _source_dir])
endif
