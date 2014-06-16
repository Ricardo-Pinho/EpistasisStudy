PROJECT (TEAM)

cmake_minimum_required(VERSION 2.6)

find_package(Boost 1.35)

link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )


SET(CMAKE_EXE_LINKER_FLAGS "-pthread")
SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS "-O3")

SET(SRC_LIST mst.cpp get_args.cpp read_data.cpp ../common/read_data.cpp test_scores.cpp parallel.cpp)

##ADD_EXECUTABLE(test test.cpp ${SRC_LIST})
##ADD_EXECUTABLE(test_origin test_origin.cpp ${SRC_LIST})

ADD_EXECUTABLE(get_snps get_snps.cpp ${SRC_LIST})
ADD_EXECUTABLE(test_all test_all.cpp ${SRC_LIST})

##ADD_EXECUTABLE(collector collector.cpp)

#ADD_EXECUTABLE(collector_all collector_all.cpp)
#ADD_EXECUTABLE(pvalue pvalue.cpp)
##ADD_EXECUTABLE(pvaluedist pvaluedist.cpp)
#ADD_EXECUTABLE(hist hist.cpp)

#target_link_libraries(./test ${Boost_LIBRARY})
#target_link_libraries(./test_orgin ${Boost_LIBRARY})

target_link_libraries(./get_snps ${Boost_LIBRARY})
target_link_libraries(./test_all ${Boost_LIBRARY})
