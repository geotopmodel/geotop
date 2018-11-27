omp = dependency('openmp', required: false)

if not omp.found()
    error('\n\nYour compiler does not support openMP. Please reconfigure with -DWITH_OMP=false\n')
endif

deps += omp
