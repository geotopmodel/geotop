
_vec =  get_option('WITH_VECTORIZATION')

if _vec == true
add_global_arguments('-march=native',
	             language : ['cpp','c'])

add_global_link_arguments('-march=native',
                          language : ['cpp','c'])

endif
              
