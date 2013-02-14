file(GLOB_RECURSE NEWTON_H "*.h")
file(GLOB_RECURSE NEWTON_CPP "*.cpp")

file(GLOB_RECURSE MAIN_CPP "main.cpp")
list(REMOVE_ITEM NEWTON_CPP ${MAIN_CPP})

if(${OLD_FREESEMIRING})
  file(GLOB_RECURSE FREESEMIRING_H "free-semiring.h")
  file(GLOB_RECURSE FREESEMIRING_CPP "free-semiring.cpp")
  file(GLOB_RECURSE FREESTRUCTURE_H "free-structure.h")
  file(GLOB_RECURSE FREESTRUCTURE_CPP "free-structure.cpp")
  list(REMOVE_ITEM NEWTON_CPP ${FREESEMIRING_H} ${FREESEMIRING_CPP}
                              ${FREESTRUCTURE_H} ${FREESTRUCTURE_CPP})
else(${OLD_FREESEMIRING})
  file(GLOB_RECURSE FREESEMIRING_H "free-semiring-old.h")
  file(GLOB_RECURSE FREESEMIRING_CPP "free-semiring-old.cpp")
  list(REMOVE_ITEM NEWTON_CPP ${FREESEMIRING_H} ${FREESEMIRING_CPP})
endif(${OLD_FREESEMIRING})

add_library(NewtonLib ${NEWTON_H} ${NEWTON_CPP})
add_executable(${PROJECTNAME} main.cpp)
target_link_libraries(${PROJECTNAME} NewtonLib ${Boost_LIBRARIES} ${LPSOLVE_LIBRARIES})
