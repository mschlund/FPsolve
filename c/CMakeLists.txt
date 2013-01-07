CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

#SET(CMAKE_C_COMPILER "clang")
#SET(CMAKE_CXX_COMPILER "clang++")

SET(PROJECTNAME newton)

PROJECT( ${PROJECTNAME} )

#FIND_PACKAGE(Boost 1.49.0 COMPONENTS program_options)
FIND_PACKAGE(Boost COMPONENTS program_options)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

FILE(GLOB_RECURSE cppsource "${PROJECT_SOURCE_DIR}/src/*.cpp")
FILE(GLOB_RECURSE csource "${PROJECT_SOURCE_DIR}/src/*.c")
FILE(GLOB_RECURSE hheaders "${PROJECT_SOURCE_DIR}/src/*.h")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")				## Optimize
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")				## Debug information
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic-errors")	## More warnings during compile
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")			## c++0x support

ADD_EXECUTABLE( ${PROJECTNAME} ${csource} ${cppsource} ${hheaders} )

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECTNAME} ${Boost_LIBRARIES}	)

