cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

project (gbons2)
find_package(CGAL REQUIRED Core)

add_subdirectory(vendor/libdts2 libdts2)

include_directories(
	${LIBDTS2_INCLUDE_DIR}
)

set(MY_C_FLAGS "${CMAKE_C_FLAGS} -frounding-math -std=c99 -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4")
set(MY_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math -std=gnu++0x -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4 -Wshadow")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")

SET(MY_LIB_LINK_LIBS
	${LIBRT_LIBRARIES}
	${LIBDTS2_LIBRARIES}
	${LIBDTS2_LINK_LIBRARIES}
	${CGAL_LIBRARIES}
	${CGAL_3RD_PARTY_LIBRARIES}
)

set(GBS_ON_S2_SOURCES_CPP
	src/gbons2.cpp
)

add_executable(${PROJECT_NAME} ${GBS_ON_S2_SOURCES_CPP})
target_link_libraries(${PROJECT_NAME} ${MY_LIB_LINK_LIBS})
add_dependencies(${PROJECT_NAME} dts2)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${LIBDTS2_COMPILE_DEFINITIONS})

