set(header_path ${CMAKE_SOURCE_DIR}/include/crsGA)

set(sources
    Logger.cpp
)

SET(headers
    ${header_path}/Export.hpp
    ${header_path}/Logger.hpp
    ${header_path}/ThreadSafeQueue.hpp
    ${header_path}/ThreadPool.hpp
    ${header_path}/Common.hpp
    ${header_path}/ScopedThread.hpp
    ${header_path}/GeneticAlgorithm.hpp
    ${header_path}/SimulatedAnnealing.hpp
)

add_library(GeneticAlgorithm STATIC ${sources} ${headers})
target_include_directories(GeneticAlgorithm
                           PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
set_target_properties(GeneticAlgorithm
                      PROPERTIES PUBLIC_HEADER "${headers}")
target_link_libraries(GeneticAlgorithm
                      PUBLIC ${CMAKE_THREAD_LIBS_INIT})


