#!/bin/bash

# # 3D tests
# rm -rf 3D/*/output-tabs-SE27XX/ 3D/*/output-maps-SE27XX/ 
# rm -rf 3D/*/geotop.log-SE27XX

# # (1) Borden05m
# ../../geotop_2.0/bin/geotop-2.0.0 3D/Borden05m 1>>3D/Borden05m/stdout.SE27XX 2>>3D/Borden05m/stderr.SE27XX
# mv 3D/Borden05m/geotop.log 3D/Borden05m/geotop.log-SE27XX
# mv 3D/Borden05m/output-tabs 3D/Borden05m/output-tabs-SE27XX
# mv 3D/Borden05m/output-maps 3D/Borden05m/output-maps-SE27XX
# mkdir 3D/Borden05m/output-tabs 3D/Borden05m/output-maps

# # (2) example
# ../../geotop_2.0/bin/geotop-2.0.0 3D/example 1>>3D/example/stdout.SE27XX 2>>3D/example/stderr.SE27XX
# mv 3D/example/geotop.log 3D/example/geotop.log-SE27XX
# mv 3D/example/output-tabs 3D/example/output-tabs-SE27XX
# mv 3D/example/output-maps 3D/example/output-maps-SE27XX
# mkdir 3D/example/output-tabs 3D/example/output-maps

# # (3) hillslope01
# ../../geotop_2.0/bin/geotop-2.0.0 3D/hillslope01 1>>3D/hillslope01/stdout.SE27XX 2>>3D/hillslope01/stderr.SE27XX
# mv 3D/hillslope01/geotop.log 3D/hillslope01/geotop.log-SE27XX
# mv 3D/hillslope01/output-tabs 3D/hillslope01/output-tabs-SE27XX
# mv 3D/hillslope01/output-maps 3D/hillslope01/output-maps-SE27XX
# mkdir 3D/hillslope01/output-tabs  3D/hillslope01/output-maps

# # (4) hillslope02_superslab
# ../../geotop_2.0/bin/geotop-2.0.0 3D/hillslope02_superslab 1>>3D/hillslope02_superslab/stdout.SE27XX 2>>3D/hillslope02_superslab/stderr.SE27XX
# mv 3D/hillslope02_superslab/geotop.log 3D/hillslope02_superslab/geotop.log-SE27XX
# mv 3D/hillslope02_superslab/output-tabs 3D/hillslope02_superslab/output-tabs-SE27XX
# mv 3D/hillslope02_superslab/output-maps 3D/hillslope02_superslab/output-maps-SE27XX
# mkdir 3D/hillslope02_superslab/output-tabs 3D/hillslope02_superslab/output-maps

# # (5) Mazia
# ../../geotop_2.0/bin/geotop-2.0.0 3D/Mazia 1>>3D/Mazia/stdout.SE27XX 2>>3D/Mazia/stderr.SE27XX
# mv 3D/Mazia/geotop.log 3D/Mazia/geotop.log-SE27XX
# mv 3D/Mazia/output-tabs 3D/Mazia/output-tabs-SE27XX
# mv 3D/Mazia/output-maps 3D/Mazia/output-maps-SE27XX
# mkdir 3D/Mazia/output-tabs 3D/Mazia/output-maps

# # (6) Muntatschini_ref_005
# ../../geotop_2.0/bin/geotop-2.0.0 3D/Muntatschini_ref_005 1>>3D/Muntatschini_ref_005/stdout.SE27XX 2>>3D/Muntatschini_ref_005/stderr.SE27XX
# mv 3D/Muntatschini_ref_005/geotop.log 3D/Muntatschini_ref_005/geotop.log-SE27XX
# mv 3D/Muntatschini_ref_005/output-tabs 3D/Muntatschini_ref_005/output-tabs-SE27XX
# mv 3D/Muntatschini_ref_005/output-maps 3D/Muntatschini_ref_005/output-maps-SE27XX
# mkdir 3D/Muntatschini_ref_005/output-tabs  3D/Muntatschini_ref_005/output-maps

# # (7) no_reflection
# ../../geotop_2.0/bin/geotop-2.0.0 3D/no_reflection 1>>3D/no_reflection/stdout.SE27XX 2>>3D/no_reflection/stderr.SE27XX
# mv 3D/no_reflection/geotop.log 3D/no_reflection/geotop.log-SE27XX
# mv 3D/no_reflection/output-tabs 3D/no_reflection/output-tabs-SE27XX
# mv 3D/no_reflection/output-maps 3D/no_reflection/output-maps-SE27XX
# mkdir 3D/no_reflection/output-tabs  3D/no_reflection/output-maps

# # (8) onepoint_hydrostatic
# ../../geotop_2.0/bin/geotop-2.0.0 3D/onepoint_hydrostatic 1>>3D/onepoint_hydrostatic/stdout.SE27XX 2>>3D/onepoint_hydrostatic/stderr.SE27XX
# mv 3D/onepoint_hydrostatic/geotop.log 3D/onepoint_hydrostatic/geotop.log-SE27XX
# mv 3D/onepoint_hydrostatic/output-tabs 3D/onepoint_hydrostatic/output-tabs-SE27XX
# mv 3D/onepoint_hydrostatic/output-maps 3D/onepoint_hydrostatic/output-maps-SE27XX
# mkdir 3D/onepoint_hydrostatic/output-tabs  3D/onepoint_hydrostatic/output-maps

# # (9) panola
# ../../geotop_2.0/bin/geotop-2.0.0 3D/panola 1>>3D/panola/stdout.SE27XX 2>>3D/panola/stderr.SE27XX
# mv 3D/panola/geotop.log 3D/panola/geotop.log-SE27XX
# mv 3D/panola/output-tabs 3D/panola/output-tabs-SE27XX
# mv 3D/panola/output-maps 3D/panola/output-maps-SE27XX
# mkdir 3D/panola/output-tabs 3D/panola/output-maps

# # (10) panola_25pixel
# ../../geotop_2.0/bin/geotop-2.0.0 3D/panola_25pixel 1>>3D/panola_25pixel/stdout.SE27XX 2>>3D/panola_25pixel/stderr.SE27XX
# mv 3D/panola_25pixel/geotop.log 3D/panola_25pixel/geotop.log-SE27XX
# mv 3D/panola_25pixel/output-tabs 3D/panola_25pixel/output-tabs-SE27XX
# mv 3D/panola_25pixel/output-maps 3D/panola_25pixel/output-maps-SE27XX
# mkdir 3D/panola_25pixel/output-tabs 3D/panola_25pixel/output-maps

# # (11) panola_25pixel_nobed
# ../../geotop_2.0/bin/geotop-2.0.0 3D/panola_25pixel_nobed 1>>3D/panola_25pixel_nobed/stdout.SE27XX 2>>3D/panola_25pixel_nobed/stderr.SE27XX
# mv 3D/panola_25pixel_nobed/geotop.log 3D/panola_25pixel_nobed/geotop.log-SE27XX
# mv 3D/panola_25pixel_nobed/output-tabs 3D/panola_25pixel_nobed/output-tabs-SE27XX
# mv 3D/panola_25pixel_nobed/output-maps 3D/panola_25pixel_nobed/output-maps-SE27XX
# mkdir 3D/panola_25pixel_nobed/output-tabs  3D/panola_25pixel_nobed/output-maps

# # (12) panola_25pixel_nobed_hydrostatic
# ../../geotop_2.0/bin/geotop-2.0.0 3D/panola_25pixel_nobed_hydrostatic 1>>3D/panola_25pixel_nobed_hydrostatic/stdout.SE27XX 2>>3D/panola_25pixel_nobed_hydrostatic/stderr.SE27XX
# mv 3D/panola_25pixel_nobed_hydrostatic/geotop.log 3D/panola_25pixel_nobed_hydrostatic/geotop.log-SE27XX
# mv 3D/panola_25pixel_nobed_hydrostatic/output-tabs 3D/panola_25pixel_nobed_hydrostatic/output-tabs-SE27XX
# mv 3D/panola_25pixel_nobed_hydrostatic/output-maps 3D/panola_25pixel_nobed_hydrostatic/output-maps-SE27XX
# mkdir 3D/panola_25pixel_nobed_hydrostatic/output-tabs 3D/panola_25pixel_nobed_hydrostatic/output-maps

# # (13) prealpiC
# ../../geotop_2.0/bin/geotop-2.0.0 3D/prealpiC 1>>3D/prealpiC/stdout.SE27XX 2>>3D/prealpiC/stderr.SE27XX
# mv 3D/prealpiC/geotop.log 3D/prealpiC/geotop.log-SE27XX
# mv 3D/prealpiC/output-tabs 3D/prealpiC/output-tabs-SE27XX
# mv 3D/prealpiC/output-maps 3D/prealpiC/output-maps-SE27XX
# mkdir 3D/prealpiC/output-tabs  3D/prealpiC/output-maps

# # (14) PSQL_test
# ../../geotop_2.0/bin/geotop-2.0.0 3D/PSQL_test 1>>3D/PSQL_test/stdout.SE27XX 2>>3D/PSQL_test/stderr.SE27XX
# mv 3D/PSQL_test/geotop.log 3D/PSQL_test/geotop.log-SE27XX
# mv 3D/PSQL_test/output-tabs 3D/PSQL_test/output-tabs-SE27XX
# mv 3D/PSQL_test/output-maps 3D/PSQL_test/output-maps-SE27XX
# mkdir 3D/PSQL_test/output-tabs 3D/PSQL_test/output-maps

# # (15) rendena
# ../../geotop_2.0/bin/geotop-2.0.0 3D/rendena 1>>3D/rendena/stdout.SE27XX 2>>3D/rendena/stderr.SE27XX
# mv 3D/rendena/geotop.log 3D/rendena/geotop.log-SE27XX
# mv 3D/rendena/output-tabs 3D/rendena/output-tabs-SE27XX
# mv 3D/rendena/output-maps 3D/rendena/output-maps-SE27XX
# mkdir 3D/rendena/output-tabs 3D/rendena/output-maps

# # (16) RottCatchment
# ../../geotop_2.0/bin/geotop-2.0.0 3D/RottCatchment 1>>3D/RottCatchment/stdout.SE27XX 2>>3D/RottCatchment/stderr.SE27XX
# mv 3D/RottCatchment/geotop.log 3D/RottCatchment/geotop.log-SE27XX
# mv 3D/RottCatchment/output-tabs 3D/RottCatchment/output-tabs-SE27XX
# mv 3D/RottCatchment/output-maps 3D/RottCatchment/output-maps-SE27XX
# mkdir 3D/RottCatchment/output-tabs 3D/RottCatchment/output-maps

# # (17) small_example
# ../../geotop_2.0/bin/geotop-2.0.0 3D/small_example 1>>3D/small_example/stdout.SE27XX 2>>3D/small_example/stderr.SE27XX
# mv 3D/small_example/geotop.log 3D/small_example/geotop.log-SE27XX
# mv 3D/small_example/output-tabs 3D/small_example/output-tabs-SE27XX
# mv 3D/small_example/output-maps 3D/small_example/output-maps-SE27XX
# mkdir 3D/small_example/output-tabs 3D/small_example/output-maps

# # (18) small_example-channel
# ../../geotop_2.0/bin/geotop-2.0.0 3D/small_example-channel 1>>3D/small_example-channel/stdout.SE27XX 2>>3D/small_example-channel/stderr.SE27XX
# mv 3D/small_example-channel/geotop.log 3D/small_example-channel/geotop.log-SE27XX
# mv 3D/small_example-channel/output-tabs 3D/small_example-channel/output-tabs-SE27XX
# mv 3D/small_example-channel/output-maps 3D/small_example-channel/output-maps-SE27XX
# mkdir 3D/small_example-channel/output-tabs  3D/small_example-channel/output-maps

# # (19) small_example-onlyEnergy
# ../../geotop_2.0/bin/geotop-2.0.0 3D/small_example-onlyEnergy 1>>3D/small_example-onlyEnergy/stdout.SE27XX 2>>3D/small_example-onlyEnergy/stderr.SE27XX
# mv 3D/small_example-onlyEnergy/geotop.log 3D/small_example-onlyEnergy/geotop.log-SE27XX
# mv 3D/small_example-onlyEnergy/output-tabs 3D/small_example-onlyEnergy/output-tabs-SE27XX
# mv 3D/small_example-onlyEnergy/output-maps 3D/small_example-onlyEnergy/output-maps-SE27XX
# mkdir 3D/small_example-onlyEnergy/output-tabs 3D/small_example-onlyEnergy/output-maps

# # (20) snow_dstr_SENSITIVITY
# ../../geotop_2.0/bin/geotop-2.0.0 3D/snow_dstr_SENSITIVITY 1>>3D/snow_dstr_SENSITIVITY/stdout.SE27XX 2>>3D/snow_dstr_SENSITIVITY/stderr.SE27XX
# mv 3D/snow_dstr_SENSITIVITY/geotop.log 3D/snow_dstr_SENSITIVITY/geotop.log-SE27XX
# mv 3D/snow_dstr_SENSITIVITY/output-tabs 3D/snow_dstr_SENSITIVITY/output-tabs-SE27XX
# mv 3D/snow_dstr_SENSITIVITY/output-maps 3D/snow_dstr_SENSITIVITY/output-maps-SE27XX
# mkdir 3D/snow_dstr_SENSITIVITY/output-tabs  3D/snow_dstr_SENSITIVITY/output-maps

# # (21) UpperAmmerCatchment
# ../../geotop_2.0/bin/geotop-2.0.0 3D/UpperAmmerCatchment 1>>3D/UpperAmmerCatchment/stdout.SE27XX 2>>3D/UpperAmmerCatchment/stderr.SE27XX
# mv 3D/UpperAmmerCatchment/geotop.log 3D/UpperAmmerCatchment/geotop.log-SE27XX
# mv 3D/UpperAmmerCatchment/output-tabs 3D/UpperAmmerCatchment/output-tabs-SE27XX
# mv 3D/UpperAmmerCatchment/output-maps 3D/UpperAmmerCatchment/output-maps-SE27XX
# mkdir 3D/UpperAmmerCatchment/output-tabs 3D/UpperAmmerCatchment/output-maps

# # (22) Vshape
# ../../geotop_2.0/bin/geotop-2.0.0 3D/Vshape 1>>3D/Vshape/stdout.SE27XX 2>>3D/Vshape/stderr.SE27XX
# mv 3D/Vshape/geotop.log 3D/Vshape/geotop.log-SE27XX
# mv 3D/Vshape/output-tabs 3D/Vshape/output-tabs-SE27XX
# mv 3D/Vshape/output-maps 3D/Vshape/output-maps-SE27XX
# mkdir 3D/Vshape/output-tabs 3D/Vshape/output-maps

# # (23) rendena
# ../../geotop_2.0/bin/geotop-2.0.0 3D/WG1_2.0_001 1>>3D/WG1_2.0_001/stdout.SE27XX 2>>3D/WG1_2.0_001/stderr.SE27XX
# mv 3D/WG1_2.0_001/geotop.log 3D/WG1_2.0_001/geotop.log-SE27XX
# mv 3D/WG1_2.0_001/output-tabs 3D/WG1_2.0_001/output-tabs-SE27XX
# mv 3D/WG1_2.0_001/output-maps 3D/WG1_2.0_001/output-maps-SE27XX
# mkdir 3D/WG1_2.0_001/output-tabs 3D/WG1_2.0_001/output-maps

# echo  3D/*/output-tabs/ | xargs -n 1 cp .placeholder
# echo  3D/*/output-tabs-SE27XX/ | xargs -n 1 cp .placeholder
# echo  3D/*/output-maps/ | xargs -n 1 cp .placeholder
# echo  3D/*/output-maps-SE27XX/ | xargs -n 1 cp .placeholder

# ls -l 3D/ | awk '{print $9}'
