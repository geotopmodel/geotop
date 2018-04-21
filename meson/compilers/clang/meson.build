_bt = get_option('buildtype') 
if _bt == 'debugoptimized' or _bt == 'release'
  add_global_arguments(['-ffp-contract=fast',
                        '-Wno-narrowing',
                        '-Wno-writable-strings'],
                       language : ['cpp','c'])
endif
