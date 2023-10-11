#!/usr/bin/bash

# get-or-build-libtorch.sh
# 
# Fetch an appropriate `libtorch`, and if necessary,
# compile it for use on this system.

# If unset, use `sudo apt-get` as the apt command
APT=${APT:-'sudo apt-get'}

# If unset, clone `pytorch` into the current directory
PYTORCH_PREFIX_PATH=${PYTORCH_PREFIX_PATH:-'.'}

# If unset, assume it was installed into the current directory
CMAKE=${CMAKE:-"$(realpath ./cmake/bin/cmake)"}

ARCH=$(uname -m)

if [ ${ARCH} == 'x86_64' ]
then

    # Download prebuilt binary
    # https://pytorch.org
    wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -q -O /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
    rm -rf /tmp/libtorch
    unzip /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -d /tmp/
    rm /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
    rm -rf ${PYTORCH_PREFIX_PATH}/pytorch
    mkdir -p ${PYTORCH_PREFIX_PATH}
    mv /tmp/libtorch ${PYTORCH_PREFIX_PATH}/pytorch

else

    # Build libtorch from sources
    cd ${PYTORCH_PREFIX_PATH}
    if [ ! -d pytorch ]; then git clone --recursive https://github.com/pytorch/pytorch; fi
    cd pytorch

    ${APT} install -y python3.10 python-is-python3 python3-pip

    _GLIBCXX_USE_CXX11_ABI=1 ${CMAKE} -S . -B build/ -DBUILD_CAFFE2=1 -DUSE_CUDA=0 -DBUILD_TEST=0 -DBUILD_PYTHON=0 -DPYTHON_EXECUTABLE=$(which python)
    make -C build/

fi

