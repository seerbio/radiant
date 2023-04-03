CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
PROJECT(PythiaDIACpp)

#INSTALL INSTRUCTIONS
# BOOST LINUX sudo apt-get install libboost-all-dev

# https://hub.docker.com/r/arrowdev/amd64-debian-10-cpp
# docker pull arrowdev/amd64-debian-10-cpp
# sudo docker run -it --entrypoint /bin/bash arrowdev/amd64-debian-10-cpp
# git clone https://github.com/apache/arrow.git
# cd arrow/cpp
# mkdir build
# cd build
# cmake .. -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DARROW_WITH_ZLIB=ON -DPARQUET_BUILD_STATIC=ON -DARROW_PARQUET=ON -DARROW_WITH_SNAPPY=ON -DPARQUET_BUILD_EXECUTABLES=ON
# make

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
SET(DLIB_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/dlib-19.24")
SET(EIGEN_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/eigen-3.4.0")
SET(ARROW_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/src")
SET(ARROW_STATIC_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libarrow.a")
SET(ARROW_STATIC_LIBRARY_DEPS "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libarrow_bundled_dependencies.a")
SET(PARQUET_STATIC_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/arrow_parquet/release/libparquet.a")
SET(ZLIB_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/zlib-1.2.12/include")
SET(ZLIB_STATIC_LIBRARY "${PROJECT_SOURCE_DIR}/ThirdPartyLibs/zlib-1.2.12/lib/libz.so")


SET(ARROW_SHARED_LIBRARY "/usr/lib/x86_64-linux-gnu/libarrow.so.1100")
SET(PARQUET_SHARED_LIBRARY "/usr/lib/x86_64-linux-gnu/libparquet.so.1100")

ADD_SUBDIRECTORY("AlgorithmsLib")
ADD_SUBDIRECTORY("Apps/IstrosLibraryBuilder")
ADD_SUBDIRECTORY("ChemLib")
ADD_SUBDIRECTORY("EigenLib")
ADD_SUBDIRECTORY("FileReadersLib")
ADD_SUBDIRECTORY("MachineLrnLib")
ADD_SUBDIRECTORY("WorkFlowsLib")
ADD_SUBDIRECTORY("UtilsLib")
#ADD_SUBDIRECTORY("ParquetLib")
