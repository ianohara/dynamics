find_package(PythonInterp REQUIRED)

set(test_dir "${CMAKE_SOURCE_DIR}/tests")
set(clar_dir "${test_dir}/clar")
set(src_dir "${CMAKE_SOURCE_DIR}/src")
set(test_executable "tests")

if (NOT EXISTS "${clar_dir}/")
    message(FATAL_ERROR "A clone of the clar repo needs to be in 'tests/'")
endif()

include_directories("${test_dir}" "${clar_dir}")

# TODO: Figure out how to fix clar so this is not necessary.  The only use of fs_copy
#       is inside an ifdef CLAR_FIXTURE_PATH.  We don't use CLAR_FIXTURE_PATH, so fs_copy
#       is defined but unused --> Boom.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function")

add_custom_command(
    OUTPUT ${test_dir}/clar.suite
    COMMAND ${PYTHON_EXECUTABLE} ${clar_dir}/generate.py -f -x${test_dir}/test.c -x${clar_dir}/ ${test_dir}
)

set_source_files_properties(
    ${clar_dir}/clar.c
    PROPERTIES OBJECT_DEPENDS ${test_dir}/clar.suite
)

add_library(clar ${clar_dir}/clar.c)

# The final test executable
add_executable(${test_executable} test.c)
target_link_libraries(${test_executable} clar)

# Add a new test.  The testsrc needs to be the .c file containing
# the tests.  The file it tests should exist in the src/ tree
# under the same name.  In other words, the tests/ directory
# should exactly mirror the src/ directory.
function(add_test testname testsrc)
    string(REPLACE " " ";" othersrc "${ARGN}")
    add_library(${testname} ${testsrc} "${src_dir}/${testsrc}" ${othersrc})
    target_link_libraries(${testname} clar c)
    target_link_libraries(${test_executable} ${testname})
endfunction()

# Add each individual test here.  The structure of the tests directory
# should exactly match the src/ directory except in the tests directory
# each .c file has tests for the corresponding src file.
add_test(test_vector "data_structures/vector.c")

add_custom_target(
    run-tests
    ALL
    COMMAND ${test_executable}
    DEPENDS ${test_executable}
)
