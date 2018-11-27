omp = dependency('openmp', required: false)

if not omp.found()
    # it might be your meson is a bit old
    compiler = meson.get_compiler('cpp')
    _c = compiler.has_argument('-fopenmp')
    _l = compiler.has_link_argument('-fopenmp')
    _h = compiler.has_header('omp.h')
    if not (_c and _l and _h)
        error('\n\nYour compiler does not support openMP. Please reconfigure with -DWITH_OMP=false\n')
    else
        omp = declare_dependency(compile_args:'-fopenmp',
                                 link_args:'-fopenmp')
    endif
endif

deps += omp
