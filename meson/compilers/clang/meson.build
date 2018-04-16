_bt = get_option('buildtype') 
if _bt == 'debugoptimized' or _bt == 'release'
  add_global_arguments('-ffp-contract=fast',
                       language : ['cpp','c'])
endif
