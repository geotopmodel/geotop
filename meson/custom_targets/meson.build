# define indent target for indenting all source files of geotop

_source_dir = meson.source_root()
_build_dir  = meson.build_root()

###########################
#    indent (astyle)      #
###########################

_indent = find_program([_source_dir+'/utils/meson_indent'], required:
                       false )

run_target('indent', command: _indent)

###########################
#    clang-tidy           #
###########################

# check that we have clang-tidy in our path
_clang_tidy = find_program('clang-tidy', required: false)

if _clang_tidy.found()
    _run_clang_tidy = find_program(_source_dir+'/utils/run_clang_tidy',
                                   required: false)

    # it is recommend to run this command from a dedicated
    # build directory defined as follows
    #
    # CXX=clang++ meson -Dbuildtype=plain build_tidy
    #
    # from build_tidy you can run
    #
    # ninja tidy
    #
    # if you want to apply the fixes suggested by the linter
    # you can do as follows
    #
    # meson configure -Dtidy_args="-fix-errors"
    # ninja tidy
    
    run_target('tidy',command: [_run_clang_tidy, _build_dir,
                                _source_dir,
                                get_option('tidy_args')])
endif
