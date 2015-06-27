file(GLOB NEWTON_TEST_H "*.h")
file(GLOB NEWTON_TEST_CPP "*.cpp")

include_directories (${PROJECT_SOURCE_DIR}/src)
link_directories (${PROJECT_BINARY_DIR}/src)

if(NOT USE_GENEPI)
	file(GLOB_RECURSE SEMILIN_SET_NDD_H "test-semilinSetNdd.h")
	file(GLOB_RECURSE SEMILIN_SET_NDD_CPP "test-semilinSetNdd.cpp")
	list(REMOVE_ITEM NEWTON_TEST_H ${SEMILIN_SET_NDD_H})
	list(REMOVE_ITEM NEWTON_TEST_CPP ${SEMILIN_SET_NDD_CPP})
endif(NOT USE_GENEPI)

if(NOT USE_LIBFA)
	file(GLOB_RECURSE LOSSY_H "test-lossy.h")
	file(GLOB_RECURSE LOSSY_CPP "test-lossy.cpp")
	list(REMOVE_ITEM NEWTON_TEST_H ${LOSSY_H})
	list(REMOVE_ITEM NEWTON_TEST_CPP ${LOSSY_CPP})
endif(NOT USE_LIBFA)

add_executable(fpsolve_test ${NEWTON_TEST_H} ${NEWTON_TEST_CPP})
target_link_libraries(fpsolve_test NewtonLib ${Boost_LIBRARIES} ${CPPUNIT_LIBRARIES} ${GENEPI_LIBRARIES} ${LIBFA_LIBRARIES} ${GMP_LIBRARIES})
