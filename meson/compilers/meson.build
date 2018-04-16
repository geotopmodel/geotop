# to do: other options


_bt = get_option('buildtype') 
if _bt == 'debugoptimized' or _bt == 'release'
  add_global_arguments(['-funroll-loops',
                        '-fstrict-aliasing',
                        '-O2'],
                       language : ['cpp','c'])
endif

_supported_compilers = ['clang', 'gcc', 'intel']

foreach _c : _supported_compilers
    if _cxx.get_id() == _c
        subdir(_c)
    endif
endforeach

