#!/usr/bin/bash

AWSSDK_PREFIX_PATH=${AWSSDK_PREFIX_PATH:-'.'}

git clone https://github.com/aws/aws-sdk-cpp.git --recurse-

#cd "$AWSSDK_PREFIX_PATH/aws-sdk-cpp" || exit
#mkdir build
#cd build || exit
#cmake .. -DBUILD_ONLY="s3" -DCMAKE_INSTALL_PREFIX=../local-install
#cmake --build . --target install

cd "$AWSSDK_PREFIX_PATH/aws-sdk-cpp" || exit
mkdir build
cd build || exit
cmake .. -DBUILD_ONLY="s3;core"
make
make install
