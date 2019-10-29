project(Banks4)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${LZ4_INCLUDE_DIR})
include_directories(${INCLUDES})

add_library(Banks4 
		SHARED 
		BBand.cpp BCalorimeter.cpp BConfig.cpp BEvent.cpp BParticle.cpp BScaler.cpp BScintillator.cpp )
target_link_libraries(Banks4 ${LZ4_LIBRARY})

#install(TARGETS Hipo4 
#		LIBRARY 
#		DESTINATION "${CMAKE_SOURCE_DIR}/lib" )
