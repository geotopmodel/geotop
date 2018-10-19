_bt = get_option('buildtype') 
if _bt == 'debugoptimized' or _bt == 'release'
  add_global_arguments('-fp_speculative=safe',
                       language : ['cpp','c'])
endif
