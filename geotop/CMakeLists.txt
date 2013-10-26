PROJECT(geotop)

	IF(COMMAND CMAKE_POLICY)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)



	ADD_EXECUTABLE(geotop
			geotop.cc
		      )

	ADD_LIBRARY(geotopPartSnow ${LINKING_MODE}
			part_snow.h
			part_snow.cc
		   )

	ADD_LIBRARY(geotopMeteo ${LINKING_MODE}
			meteo.h
			meteo.cc
			meteodistr.cc
			meteodistr.h
			meteodata.cc
			meteodata.h
			times.cc
			times.h
		   )

	ADD_LIBRARY(geotopCommon ${LINKING_MODE}
			constants.h
			datastructs.h
			geotop_common.cc
			geotop_common.h
			keywords.h
			keywords.cc
		   )

	ADD_LIBRARY(geotopEnergyBalance ${LINKING_MODE}
			energy.balance.cc
			energy.balance.h
			turbulence.cc
			turbulence.h
			snow.cc
			snow.h
			radiation.cc
			radiation.h
			pedo.funct.cc
			pedo.funct.h
			vegetation.cc
			vegetation.h
			clouds.cc
			clouds.h
		   )

	ADD_LIBRARY(geotopL ${LINKING_MODE}

			PBSM.cc
			PBSM.h
			blowingsnow.cc
			blowingsnow.h
			channels.cc
			channels.h
			deallocate.cc
			deallocate.h
			energy.balance.h
			indices.cc
			indices.h
			input.cc
			input.h
			output.cc
			output.h
			output_nc.cc
	output_nc.h
	parameters.cc
	parameters.h
	pure_gt.cc
	pure_gt.h
	recovering.cc
	recovering.h
	struct.geotop.h
	tables.cc
	tables.h
	water.balance.cc
	water.balance.h
	)

	TARGET_LINK_LIBRARIES(geotopCommon
			${FLUIDTURTLE_LIBRARIES}
#${ASCII_LIBRARIES}
#${METEOIO_PLUGIN_INCLUDE_PATH}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_PART_SNOW_LIBRARIES}
			${GEOTOP_COMMON_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_METEO_LIBRARIES}
			${GEOTOP_COMMON_LIBRARIES}
			${ASCII_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_ENERGY_BALANCE_LIBRARIES}
			${GEOTOP_METEO_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_LIBRARIES}
			${GEOTOP_ENERGY_BALANCE_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(geotop
			${GEOTOP_LIBRARIES}
			)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

