file(GLOB_RECURSE NEWTON_H "*.h")
file(GLOB_RECURSE NEWTON_CPP "*.cpp")

file(GLOB_RECURSE MAIN_CPP "main.cpp")
list(REMOVE_ITEM NEWTON_CPP ${MAIN_CPP})

if(${OLD_SEMILINEAR_SET})
  file(GLOB_RECURSE SPARSE_VEC_H "sparse_vec.h")
  file(GLOB_RECURSE LINEAR_SET_H "linear_set.h")
  file(GLOB_RECURSE SEMILINEAR_SET_H "semilinear_set.h")
  list(REMOVE_ITEM NEWTON_H ${SPARSE_VEC_H} ${LINEAR_SET_H} ${SEMILINEAR_SET_H})
else(${OLD_SEMILINEAR_SET})
  file(GLOB_RECURSE SEMILIN_SET_EXP_H "semilinSetExp.h")
  file(GLOB_RECURSE SEMILIN_SET_EXP_CPP "semilinSetExp.cpp")
  list(REMOVE_ITEM NEWTON_H ${SEMILIN_SET_EXP_H})
  list(REMOVE_ITEM NEWTON_CPP ${SEMILIN_SET_EXP_CPP})
endif(${OLD_SEMILINEAR_SET})

set_source_files_properties(parser.cpp PROPERTIES COMPILE_FLAGS -g0)

add_library(NewtonLib ${NEWTON_H} ${NEWTON_CPP})
add_executable(${PROJECTNAME} main.cpp)
target_link_libraries(${PROJECTNAME} NewtonLib ${Boost_LIBRARIES} ${LPSOLVE_LIBRARIES} fa rt)