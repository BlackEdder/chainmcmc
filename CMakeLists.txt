cmake_minimum_required(VERSION 2.6)
project(chainmcmc)

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall")
#SET (CMAKE_CXX_FLAGS "-g -std=c++11 -Wall")
#SET (CMAKE_CXX_FLAGS "-O2 -pipe -ffast-math -fomit-frame-pointer -Wfatal-errors -std=gnu++0x")

#link_directories( "/usr/lib/R/lib" )
find_library( CPPA cppa )
if (NOT CPPA)
	message( FATAL_ERROR "libcppa not found" )
	endif()

include_directories("include")

SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
add_library( chainmcmc SHARED
	src/chainmcmc/chain.cc;src/chainmcmc/trace.cc;src/chainmcmc/file.cc;src/chainmcmc/logger.cc;src/chainmcmc/prior.cc;src/chainmcmc/temperature.cc)
target_link_libraries( chainmcmc ${CPPA})

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
FILE(GLOB examplesRaw RELATIVE ${CMAKE_HOME_DIRECTORY} examples/*.cc)
FOREACH(example ${examplesRaw})
	get_filename_component( ex ${example} NAME_WE )
	add_executable(${ex} EXCLUDE_FROM_ALL ${example})
 	target_link_libraries(${ex} chainmcmc)
 	SET(examples ${examples};${ex})
ENDFOREACH()

add_custom_target( examples
	DEPENDS ${examples}
	)


#Install the library
install (TARGETS chainmcmc DESTINATION lib)
install (FILES include/chainmcmc/chainmcmc.hh
	include/chainmcmc/chain.hh
	include/chainmcmc/file.hh
	include/chainmcmc/trace.hh
	include/chainmcmc/logger.hh
	include/chainmcmc/prior.hh
	include/chainmcmc/parameter.hh
	include/chainmcmc/temperature.hh
 	DESTINATION include/chainmcmc)

#Documentation
set(DOXYFILE_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/api-doc/chainmcmc")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(UseDoxygen)


#Test
if (CMAKE_BUILD_TYPE MATCHES Debug)
find_package(CxxTest)
if(CXXTEST_FOUND)
	enable_testing()
	FILE(GLOB testsRaw RELATIVE ${CMAKE_HOME_DIRECTORY} tests/test_*.h)
	FOREACH(test ${testsRaw})
		get_filename_component( ftest ${test} NAME_WE )
		CXXTEST_ADD_TEST(${ftest} ${ftest}.cc ${test})
		target_link_libraries(${ftest} chainmcmc)
	ENDFOREACH()
endif()
endif()

