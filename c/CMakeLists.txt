CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

SET(PROJECTNAME newton)

PROJECT( ${PROJECTNAME} )

FIND_PACKAGE(Boost 1.49.0 COMPONENTS program_options)

IF(EXISTS "src/")
        FILE(GLOB_RECURSE cppsource "src/*.cpp")
        FILE(GLOB_RECURSE csource "src/*.c")
        FILE(GLOB_RECURSE hheaders "src/*.h")
ELSE()
        FILE(GLOB cppsource "*.cpp")
        FILE(GLOB csource "*.c")
        FILE(GLOB hheaders "*.h")
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")        ## Optimize
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")         ## Debug information
        #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")      ## More warnings during compile
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x") ## c++0x support
ENDIF()

ADD_EXECUTABLE( ${PROJECTNAME} ${csource} ${cppsource} ${hheaders} )

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECTNAME} ${Boost_LIBRARIES}	)

