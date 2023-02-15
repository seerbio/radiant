#cmake_minimum_required(VERSION 3.24)
#project(PythiaDIACpp)
#
#set(CMAKE_CXX_STANDARD 17)
#
#set(SOURCE_FILES
#        main.cpp
#        src/DeisotoperTandem.cpp
#        src/EigenSparseUtils.cpp
#        src/Error.cpp
#        src/ErrorUtils.cpp
#        src/GlobalSettings.cpp
#        src/MathUtils.cpp
#        src/ParquetReader.cpp
#        src/StringUtils.cpp
#       )
#
#set(HEADER_FILES
#        src/ClassAliases.h
#        src/DeisotoperTandem.h
#        src/EigenSparseUtils.h
#        src/Error.h
#        src/ErrorUtils.h
#        src/GlobalSettings.h
#        src/MathUtils.h
#        src/ParquetReader.h
#        src/StringUtils.h
#        )
#
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#
#set(EIGEN_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/eigen-3.4.0")
#message(Eigen Dir, ${EIGEN_DIRECTORY})
#
#set(ARROW_SHARED_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/ParquetFiles/libarrow.so")
#set(PARQUET_SHARED_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/ParquetFiles/libparquet.so")
#
#include_directories(
#        "${EIGEN_DIRECTORY}"
#)
#message(Arrow Include Directory: ${ARROW_INCLUDE_DIR})
#
#add_executable(SparkDIA main.cpp ${SOURCE_FILES} ${HEADER_FILES})
#target_link_libraries(SparkDIA PRIVATE ${ARROW_SHARED_LIBRARY} ${PARQUET_SHARED_LIBRARY})
