CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
PROJECT(PythiaDIACpp)

list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/pytorch/build/share/cmake/Torch")
set(CMAKE_TOOLCHAIN_FILE "${PROJECT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/aws-sdk-cpp/build")
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/aws-sdk-cpp/cmake")
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/aws-sdk-cpp/build/generated/src/aws-cpp-sdk-s3")
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/aws-sdk-cpp/build/src/aws-cpp-sdk-core")
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/aws-sdk-cpp/crt/aws-crt-cpp/cmake")
#list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/aws-sdk-cpp/crt/aws-crt-cpp/crt")

#function(recurse_subdirectories base_path)
#    file(GLOB_RECURSE subdirs RELATIVE ${base_path} ${base_path}/*)
#    foreach(subdir ${subdirs})
#        if(IS_DIRECTORY "${base_path}/${subdir}")
#            list(APPEND CMAKE_PREFIX_PATH "${base_path}/${subdir}")
#        endif()
#    endforeach()
#endfunction()
#
## Initial base path to start the search
#set(BASE_PATH "${PROJECT_SOURCE_DIR}/aws-sdk-cpp/)
#
## Populate CMAKE_PREFIX_PATH
#recurse_subdirectories(${BASE_PATH})

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

ADD_SUBDIRECTORY("AlgorithmsFFLib")
ADD_SUBDIRECTORY("Apps/LycaonMzMLTransmorgitron")
ADD_SUBDIRECTORY("Apps/PythiaDIA")
ADD_SUBDIRECTORY("AWSLib")
ADD_SUBDIRECTORY("ChemLib")
ADD_SUBDIRECTORY("EigenLib")
ADD_SUBDIRECTORY("FileReadersLib")
ADD_SUBDIRECTORY("PyTorchLib")
ADD_SUBDIRECTORY("WorkFlowsLib")
ADD_SUBDIRECTORY("UtilsLib")
