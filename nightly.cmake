##This is a script file to automate compilation and execution of
##DynamO tests. This is used for continuous integration testing.

set(BRANCH "experimental")
#Please set the following to a description of the test system (i.e. CentOS6.5)
set(CTEST_SITE "UserSite")
set(CTEST_BUILD_NAME "linux-gcc-default")

##########################################################################
# Site-specific setup
##########################################################################
set(CTEST_SOURCE_DIRECTORY "/home/mjki2mb2/CTest/dynamo/${BRANCH}/source")
set(CTEST_BINARY_DIRECTORY "/home/mjki2mb2/CTest/dynamo/${BRANCH}/build")

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_CONFIGURATION "Release")
set(CTEST_BUILD_FLAGS "-j1")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
SET (CTEST_COMMAND "ctest")
find_program(CTEST_GIT_COMMAND NAMES git)

#set(WITH_MEMCHECK TRUE)
#set(WITH_COVERAGE TRUE)
#find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
#find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/tests/valgrind.supp)

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone -b ${BRANCH} https://github.com/toastedcrumpets/DynamO.git ${CTEST_SOURCE_DIRECTORY}")
endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTING:BOOL=ON \"-DCTEST_BUILD_FLAGS:STRING=${CTEST_BUILD_FLAGS}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

while (1)
 set(START_TIME ${CTEST_ELAPSED_TIME})
 ctest_start("Continuous")
 ctest_update(RETURN_VALUE HAD_UPDATES)
 if(${HAD_UPDATES} GREATER 0)
  ctest_submit(PARTS Update Notes)
  ctest_configure()
  ctest_submit(PARTS Update Notes Configure)
  ctest_build()
  ctest_submit(PARTS Update Notes Configure Build)
  ctest_test()
  ctest_submit(PARTS Update Notes Configure Build Test)
 endif()
 ctest_sleep( ${START_TIME} 60 ${CTEST_ELAPSED_TIME})
endwhile()
