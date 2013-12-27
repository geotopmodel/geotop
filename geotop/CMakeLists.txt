PROJECT(geotop)

IF(COMMAND CMAKE_POLICY)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


	ADD_EXECUTABLE(geotop
			geotop.cc
		      )

	ADD_LIBRARY(${GEOTOP_PART_SNOW_LIBRARIES} ${LINKING_MODE}
			part_snow.h
			part_snow.cc
		   )

	ADD_LIBRARY(${GEOTOP_METEO_LIBRARIES} ${LINKING_MODE}
			meteo.h
			meteo.cc
			meteodistr.cc
			meteodistr.h
			meteodata.cc
			meteodata.h
			times.cc
			times.h

			snow.cc
			snow.h
			turbulence.cc
			turbulence.h
			radiation.cc
			radiation.h
			vegetation.cc
			vegetation.h
			clouds.cc
			clouds.h

			pedo.funct.cc
			pedo.funct.h
		   )

	ADD_LIBRARY(${INPUT_KEYWORDS_LIBRARIES} ${LINKING_MODE}
			inputKeywords.cc
			inputKeywords.h
		   )

	ADD_LIBRARY(${GEOTOP_COMMON_LIBRARIES} ${LINKING_MODE}
			constants.h
			datastructs.h
			struct.geotop.h
			geotop_common.cc
			geotop_common.h
			keywords.h
			keywords.cc
		   )

	ADD_LIBRARY(${GEOTOP_ENERGY_BALANCE_LIBRARIES} ${LINKING_MODE}
			energy.balance.cc
			energy.balance.h
		   )


	ADD_LIBRARY(${GEOTOP_WATER_BALANCE_LIBRARIES} ${LINKING_MODE}
			water.balance.cc
			water.balance.h
	)

	ADD_LIBRARY(${GEOTOP_LIBRARIES} ${LINKING_MODE}
            PBSM.cc
			PBSM.h
			blowingsnow.cc
			blowingsnow.h
			channels.cc
			channels.h
			indices.cc
			indices.h
			input.cc
			input.h
            output_nc.h
            parameters.cc
            parameters.h
            pure_gt.cc
            pure_gt.h
            recovering.cc
            recovering.h
	)

	ADD_LIBRARY(geotopOutPostProcessing ${LINKING_MODE}
			tables.cc
			tables.h
			output.cc
			output_nc.cc
			output.h
			deallocate.cc
			deallocate.h
	)

	TARGET_LINK_LIBRARIES(${GEOTOP_COMMON_LIBRARIES}
			${FLUIDTURTLE_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_PART_SNOW_LIBRARIES}
			${GEOTOP_COMMON_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_OUT_POSTPROC_LIBRARIES}
			${GEOTOP_METEO_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_METEO_LIBRARIES}
			${GEOTOP_COMMON_LIBRARIES}
			${ASCII_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_WATER_BALANCE_LIBRARIES}
			${GEOTOP_METEO_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_ENERGY_BALANCE_LIBRARIES}
			${GEOTOP_METEO_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(${INPUT_KEYWORDS_LIBRARIES}
            ${Boost_FILESYSTEM_LIBRARY}
            ${Boost_THREAD_LIBRARY}
            ${Boost_SYSTEM_LIBRARY}
            ${Boost_IOSTREAMS_LIBRARY}
            ${Boost_REGEX_LIBRARY}
			)

	TARGET_LINK_LIBRARIES(${GEOTOP_LIBRARIES}
            ${INPUT_KEYWORDS_LIBRARIES}
			${GEOTOP_OUT_POSTPROC_LIBRARIES}
			${GEOTOP_ENERGY_BALANCE_LIBRARIES}
			${GEOTOP_WATER_BALANCE_LIBRARIES}
			)

	TARGET_LINK_LIBRARIES(geotop
			${GEOTOP_LIBRARIES}
			)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

