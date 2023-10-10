CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
PROJECT(PythiaDIACpp)

#INSTALL INSTRUCTIONS
# BOOST LINUX sudo apt-get install libboost-all-dev

#PYTORCH
# instructions at: https://pytorch.org/cppdocs/installing.html
# c api: https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
# after unzipping, move the entire libtorch folder to ThirdPartyLibs folder

# ARROW AND PARQUET
# https://hub.docker.com/r/arrowdev/amd64-debian-10-cpp
# docker pull arrowdev/amd64-debian-10-cpp
# sudo docker run -it --entrypoint /bin/bash arrowdev/amd64-debian-10-cpp
# git clone https://github.com/apache/arrow.git
# cd arrow/cpp
# mkdir build
# cd build
# cmake .. -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DARROW_WITH_ZLIB=ON -DPARQUET_BUILD_STATIC=ON -DARROW_PARQUET=ON -DARROW_WITH_SNAPPY=ON -DPARQUET_BUILD_EXECUTABLES=ON
# make

### INSTRUCTIONS TO PULL FILES FROM DOCKER CONTAINER
# zip the release directory

# Open separate termial window
# sudo docker ps
# sudo docker exec <container_id> chmod 777 /path/to/zipped/directory

# docker cp <container_name>:path/to/file -
# docker cp <container_name>:path/to/file - > file.txt

# sudo docker cp 7449c6433401:/arrow/cpp/arrow.zip -
# sudo docker cp 7449c6433401:/arrow/cpp/arrow.zip - > arrow.zip

# sudo docker cp d99e58c57b48:/arrow/cpp/build/arrow_parquet.zip -
# sudo docker cp d99e58c57b48:/arrow/cpp/build/arrow_parquet.zip - > arrow_parquet.zip
########################################################

# LIBRARIES TO ADD TO FRESH INSTALL
# sudo apt-get install libutf8proc2
# sudo apt-get install libre2-5
# sudo apt-get install libsnappy1v5
# sudo apt-get install libthrift-0.13.0

enable_testing(true)

SET(ROOT_PROJECT_NAME ${PROJECT_NAME})

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTORCC ON)
SET(CMAKE_AUTOUIC ON)

SET(CMAKE_INCLUDE_CURRENT_DIR ON)

SET(FUNCTIONAL_PLUS_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/FunctionalPlus/include")
SET(FRUGALLY_DEEP_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/frugally-deep/include")
SET(JSON_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/json/include")
SET(MODELS_DIRECTORY "${PROJECT_SOURCE_DIR}/MachineLrnLib/Models")
SET(EIGEN_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/eigen-3.4.0")
SET(ARROW_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_src")
SET(ARROW_STATIC_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libarrow.a")
SET(ARROW_STATIC_LIBRARY_DEPS "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libarrow_bundled_dependencies.a")
SET(PARQUET_STATIC_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libparquet.a")
SET(ARROW_SHARED_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libarrow.so.1100")
SET(PARQUET_SHARED_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libparquet.so.1100")
#SET(PYTORCH_PREFIX_PATH "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/libtorch")
SET(PYTORCH_PREFIX_PATH "/libtorch")

ADD_SUBDIRECTORY("AlgorithmsLib")
ADD_SUBDIRECTORY("Apps/IstrosLibraryBuilder")
ADD_SUBDIRECTORY("Apps/LycaonMzMLTransmorgitron")
ADD_SUBDIRECTORY("Apps/OrpheusFastaShredder")
ADD_SUBDIRECTORY("Apps/PythiaDIA")
ADD_SUBDIRECTORY("ChemLib")
ADD_SUBDIRECTORY("EigenLib")
ADD_SUBDIRECTORY("FileReadersLib")
ADD_SUBDIRECTORY("KarnnLib")
ADD_SUBDIRECTORY("MachineLrnLib")
ADD_SUBDIRECTORY("PyTorchLib")
ADD_SUBDIRECTORY("WorkFlowsLib")
ADD_SUBDIRECTORY("UtilsLib")


## Add a custom target for running all tests
#add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND})
#
## Add a test for each QTest executable
#file(GLOB_RECURSE test_executables ${CMAKE_CURRENT_BINARY_DIR}/*Test)
#foreach(test_executable ${test_executables})
#    get_filename_component(test_name ${test_executable} NAME_WE)
#    add_test(NAME ${test_name} COMMAND ${test_executable})
#    add_dependencies(run_tests ${test_name})
#endforeach()

# /snap/clion/234/bin/cmake/linux/x64/bin/cmake --build /home/anichols/Repositories/Builds/PythiaDIACpp --target run_tests -j 6