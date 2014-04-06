file(GLOB_RECURSE NEWTON_H "*.h")
file(GLOB_RECURSE NEWTON_CPP "*.cpp")

file(GLOB_RECURSE MAIN_CPP "main.cpp")
list(REMOVE_ITEM NEWTON_CPP ${MAIN_CPP})

if(NOT USE_GENEPI)
	file(GLOB_RECURSE SEMILIN_SET_NDD_H "semilinSetNdd*.h")
	file(GLOB_RECURSE SEMILIN_SET_NDD_CPP "semilinSetNdd*.cpp")
	list(REMOVE_ITEM NEWTON_H ${SEMILIN_SET_NDD_H})
	list(REMOVE_ITEM NEWTON_CPP ${SEMILIN_SET_NDD_CPP})
endif(NOT USE_GENEPI)
if(USE_GENEPI)
	set_source_files_properties(parser.cpp main.cpp PROPERTIES COMPILE_DEFINITIONS "USE_GENEPI")
endif(USE_GENEPI)

set_source_files_properties(parser.cpp PROPERTIES COMPILE_FLAGS -g0)

add_library(NewtonLib ${NEWTON_H} ${NEWTON_CPP})
add_executable(${PROJECTNAME} main.cpp)
target_link_libraries(${PROJECTNAME} NewtonLib ${Boost_LIBRARIES} ${LPSOLVE_LIBRARIES} ${GENEPI_LIBRARIES} )
