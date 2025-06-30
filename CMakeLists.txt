CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
PROJECT(PythiaDIACpp)

enable_testing(true)

SET(ROOT_PROJECT_NAME ${PROJECT_NAME})

SET(CMAKE_CXX_STANDARD 17)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

# ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

SET(EIGEN_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/eigen-3.4.90")
SET(NANOFLANN_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/nanoflann/include")
SET(TIMS_INCLUDE_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/timsdata-2.21.0.4/timsdata/include/c/")
SET(TIMS_CPPSQL_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/timsdata-2.21.0.4/timsdata/thirdparty/CppSQLite-master")
SET(TIMS_CPPSQL_HEADER "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/timsdata-2.21.0.4/timsdata/thirdparty/CppSQLite-master/CppSQLite3.h")
SET(TIMS_CPPSQL_SOURCE "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/timsdata-2.21.0.4/timsdata/thirdparty/CppSQLite-master/CppSQLite3.cpp")
SET(TIMS_LIB "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/timsdata-2.21.0.4/timsdata/linux64/libtimsdata.so")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto=auto")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto=auto")

ADD_SUBDIRECTORY("AlgorithmsFFLib")
ADD_SUBDIRECTORY("Apps/LycaonMzMLTransmorgitron")
ADD_SUBDIRECTORY("Apps/FeatureOptimizerPrime")
ADD_SUBDIRECTORY("Apps/PythiaDIA")
ADD_SUBDIRECTORY("ChemLib")
ADD_SUBDIRECTORY("EigenLib")
ADD_SUBDIRECTORY("FileReadersLib")
ADD_SUBDIRECTORY("PyTorchLib")
ADD_SUBDIRECTORY("WorkFlowsLib")
ADD_SUBDIRECTORY("UtilsLib")
