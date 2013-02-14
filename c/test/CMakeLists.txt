file(GLOB_RECURSE NEWTON_TEST_H "*.h")
file(GLOB_RECURSE NEWTON_TEST_CPP "*.cpp")

include_directories (${PROJECT_SOURCE_DIR}/src)
link_directories (${PROJECT_BINARY_DIR}/src)

add_executable(newton_test ${NEWTON_TEST_H} ${NEWTON_TEST_CPP})
target_link_libraries(newton_test NewtonLib ${Boost_LIBRARIES} ${CPPUNIT_LIBRARIES} ${LPSOLVE_LIBRARIES})
