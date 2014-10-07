##########################################################################
# Site-specific setup
##########################################################################

# Dashboard model (Continuous, Experimental, Nightly)
set(CTEST_SOURCE_DIRECTORY "/home/mjki2mb2/CTest/dynamo/source")
set(CTEST_BINARY_DIRECTORY "/home/mjki2mb2/CTest/dynamo/build")

set(CTEST_SITE "CentOS6.5")
set(CTEST_BUILD_NAME "linux-gcc-default")

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CF_BUILD_CONFIGURATION "Release")

set(CTEST_BUILD_FLAGS "-j8")
set(CTEST_BUILD_OPTIONS "\"-DCTEST_BUILD_FLAGS:STRING=${CTEST_BUILD_FLAGS}\"")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
SET (CTEST_COMMAND "ctest")
find_program(CTEST_GIT_COMMAND NAMES git)

#set(WITH_MEMCHECK TRUE)
#set(WITH_COVERAGE TRUE)
#find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
#find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/tests/valgrind.supp)

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone https://github.com/toastedcrumpets/DynamO.git ${CTEST_SOURCE_DIRECTORY}")
endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTING:BOOL=ON ${CTEST_BUILD_OPTIONS}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

while (1)
 set(START_TIME ${CTEST_ELAPSED_TIME})
 ctest_start("Continuous")
 ctest_update(RETURN_VALUE HAD_UPDATES)
 if(${HAD_UPDATES} GREATER 0)
  ctest_configure()
  ctest_build()
  ctest_test()
  ctest_submit()
 endif()
 ctest_sleep( ${START_TIME} 60 ${CTEST_ELAPSED_TIME})
endwhile()

##########################################################################
# Configuration
##########################################################################

# Initial cache entries
set(CF_CONFIG_OPTIONS "\"-DCTEST_BUILD_FLAGS:STRING=${CTEST_BUILD_FLAGS}\" -DCMAKE_BUILD_TYPE:STRING=${CF_BUILD_TYPE}")
set(CF_CONFIG_OPTIONS "${CF_CONFIG_OPTIONS} -DSITE:STRING=${CTEST_SITE}")

# Assemble cmake command line
set(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND} ${CF_CONFIG_OPTIONS}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

##########################################################################
# Run the dashboard continuously
##########################################################################
while (True)
 set(START_TIME ${CTEST_ELAPSED_TIME})
 message(STATUS "Starting CTest")
 ctest_start(${MODEL})
 message(STATUS "Updating code")
 ctest_update(RETURN_VALUE HAD_UPDATES)
 message(STATUS "Checking for updates")
 if(${HAD_UPDATES} GREATER 0)
    message(STATUS "Updates found, beginning build")
   ctest_configure()
   ctest_build()
   ctest_test()
   ctest_submit()
 else()
    message(STATUS "No updates")
 endif()
 ctest_sleep( ${START_TIME} 60 ${CTEST_ELAPSED_TIME})
endwhile()