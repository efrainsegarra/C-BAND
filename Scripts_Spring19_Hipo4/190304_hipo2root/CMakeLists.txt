set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include )


include_directories(${INCLUDES})
add_executable(hipo2root hipo2root.cpp )
target_link_libraries(hipo2root ${ROOT_LIBRARIES} Hipo4 Banks4 )
