#!/usr/bin/bash

# get-or-build-libtorch.sh
# 
# Fetch an appropriate `libtorch`, and if necessary,
# compile it for use on this system.

PYTORCH_PREFIX_PATH=${PYTORCH_PREFIX_PATH:-'.'}

# https://pytorch.org
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -q -O /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
rm -rf /tmp/libtorch
unzip /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -d /tmp/
rm /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
rm -rf ${PYTORCH_PREFIX_PATH}/pytorch
mv /tmp/libtorch ${PYTORCH_PREFIX_PATH}/pytorch

## Build libtorch from sources
#WORKDIR /src/
#RUN git clone --recursive https://github.com/pytorch/pytorch \
#    && cd pytorch \
#    && _GLIBCXX_USE_CXX11_ABI=1 python tools/build_libtorch.py

