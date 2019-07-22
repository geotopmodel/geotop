_bt = get_option('buildtype') 
if _bt == 'debugoptimized' or _bt == 'release'
  add_global_arguments(['-ffp-contract=fast'],
                       language : ['cpp','c'])
endif

add_global_arguments(['-Wno-narrowing',
                      '-Wno-writable-strings',
                      '-Wno-c++11-narrowing'],
                     language : ['cpp','c'])
