
set(_tests
  1D/ARF_1D
  1D/B2_BeG_017
  1D/Bro 1D/Calabria
  1D/ColdelaPorte 
  1D/CostantMeteo 
  1D/Jungfraujoch 
  1D/Matsch_B2_Ref_007 
  1D/Matsch_P2_Ref_007 
  1D/PureDrainage 
  1D/PureDrainageFaked 
  1D/PureDrainageRainy 
  1D/PureDrainageRainySlope 
  3D/Borden05m 
  3D/example 
  3D/hillslope01 
  3D/hillslope02_superslab 
  3D/Mazia
  3D/Muntatschini_ref_005
  3D/no_reflection 
  3D/onepoint_hydrostatic 
  3D/panola 
  3D/panola_25pixel 
  3D/panola_25pixel_nobed 
  3D/panola_25pixel_nobed_hydrostatic 
  3D/prealpiC 
  3D/PSQL_test 
  3D/rendena
  3D/RottCatchment
  3D/small_example 
  3D/small_example-channel 
  3D/small_example-onlyEnergy 
  3D/snow_dstr_SENSITIVITY
  3D/UpperAmmerCatchment
  3D/Vshape 
  3D/WG1_2.0_001)

foreach(_t ${_tests})
        add_test(NAME "${_t}"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_t}
                COMMAND $<TARGET_FILE:geotop> .)
endforeach()

# setup the test_runner tests. They have an additional suite entry,
# namely 'test_runner'

foreach(_t ${_tests})
  add_test(NAME "${_t}.test_runner"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_t}
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/compare.sh)
  set_tests_properties("${_t}.test_runner" PROPERTIES
	DEPENDS "${_t}")
endforeach()
