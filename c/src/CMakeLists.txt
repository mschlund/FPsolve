file(GLOB_RECURSE NEWTON_H "*.h")
file(GLOB_RECURSE NEWTON_CPP "*.cpp")
file(GLOB_RECURSE MAIN_CPP "main.cpp")

list(REMOVE_ITEM NEWTON_CPP ${MAIN_CPP})

add_library(NewtonLib ${NEWTON_H} ${NEWTON_CPP})
add_executable(${PROJECTNAME} main.cpp)
target_link_libraries(${PROJECTNAME} NewtonLib ${Boost_LIBRARIES})
