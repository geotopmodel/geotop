# this variable is useful when finding libraries and for setting the
# proper flags for the compiler
_cxx = meson.get_compiler('cpp')

# find and handle all the dependencies
subdir('dependencies')

# set proper compile flags
subdir('compilers')

# generate config.h and version.h
subdir('config')

# custom targets
subdir('custom_targets')

# install target
subdir('install')
