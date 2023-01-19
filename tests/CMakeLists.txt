set(PackageName ${PROJECT_NAME})

configure_file(${DD4hep_ROOT}/cmake/run_test_package.sh run_test_${PackageName}.sh @ONLY)
install(PROGRAMS
  ${CMAKE_CURRENT_BINARY_DIR}/run_test_${PackageName}.sh
  DESTINATION bin)

# Helper function to quickly write tests for luxegeo. It takes care of setting
# up the correct environment and also to register it with CTest. The resulting
# tests will fail if the command fails or if it prints anything pointing towards
# a runtime error to screen
#
# Arguments:
#   name        The name of the test (as it will also appear in the CTest output)
#   command     The command to run
function(LUXEGEO_TEST name)
  add_test(${name} "${CMAKE_INSTALL_PREFIX}/bin/run_test_${PackageName}.sh"
    ${ARGN})
  set_tests_properties(${name}
    PROPERTIES
      FAIL_REGULAR_EXPRESSION "Exception;EXCEPTION;ERROR;Error"
    )
endfunction()

LUXEGEO_TEST(test_SimpleTracker ddsim --compactFile ${CMAKE_SOURCE_DIR}/compact/SimpleTracker.xml -G -N 1 --gun.multiplicity 100 --gun.particle e+)

LUXEGEO_TEST(test_LUXETracker ddsim --compactFile ${CMAKE_SOURCE_DIR}/compact/LUXETracker.xml -G -N 1 --gun.multiplicity 100 --gun.particle e+)
