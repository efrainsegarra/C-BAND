project(Hipo4)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${LZ4_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_library(Hipo4 
		bank.cpp dictionary.cpp event.cpp reader.cpp record.cpp recordbuilder.cpp utils.cpp wrapper.cpp writer.cpp )

target_link_libraries(Hipo4 ${LZ4_LIBRARY} )
