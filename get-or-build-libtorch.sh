#!/usr/bin/bash

# get-or-build-libtorch.sh
# 
# Fetch an appropriate `libtorch`, and if necessary,
# compile it for use on this system.
#
# If libtorch is built from sources, the repository
# will be used or created at `${PYTORCH_PREFIX_PATH}/pytorch`
# (if not specified, `${PYTORCH_PREFIX_PATH}` will be the
# working directory).
#
# If a pre-built libtorch is used, it will be placed in
# `${PYTORCH_PREFIX_PATH}/pytorch/build` to match the location
# of files when building from sources.

# If unset, use `sudo apt-get` as the apt command
APT=${APT:-'sudo apt-get'}

# If unset, place `pytorch` into the current directory
PYTORCH_PREFIX_PATH=${PYTORCH_PREFIX_PATH:-'.'}

# If unset, assume it was installed into the current directory;
# see `install-build-deps-ubuntu.sh`.
CMAKE=${CMAKE:-"$(realpath ./cmake/bin/cmake)"}

ARCH=$(uname -m)
if [ "${ARCH}" == 'x86_64' ]
then

    echo "Found x86 architecture (${ARCH})"
    echo "Using prebuilt libtorch..."

    # Download prebuilt binary
    # https://pytorch.org
    wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -q -O /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
    rm -rf /tmp/libtorch
    unzip /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -d /tmp/
    rm /tmp/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
    rm -rf "${PYTORCH_PREFIX_PATH}/pytorch/build/"
    mkdir -p "${PYTORCH_PREFIX_PATH}/pytorch/"
    mv /tmp/libtorch "${PYTORCH_PREFIX_PATH}/pytorch/build"

    echo "Done copying libtorch to ${PYTORCH_PREFIX_PATH}/pytorch/build/"

else

    echo "Found architecture: ${ARCH}"
    echo "Building libtorch..."

    # Install Python, which is needed by the libtorch build
    ${APT} install -y python3.10 python-is-python3 python3-pip git
    pip install setuptools pyyaml typing-extensions

    # Build libtorch from sources
    mkdir -p "${PYTORCH_PREFIX_PATH}"
    cd "${PYTORCH_PREFIX_PATH}" || exit
    if [ ! -d pytorch ]; then git clone --recursive https://github.com/pytorch/pytorch; fi
    cd pytorch || exit

    if [ -f "${CMAKE}" ]
    then
        CMAKE=$(realpath  "${CMAKE}")
        CMAKE_DIR=$(dirname "${CMAKE}")
        PATH="${CMAKE_DIR}:${PATH}"
    fi
    #_GLIBCXX_USE_CXX11_ABI=1 USE_CUDA=OFF USE_OPENMP=OFF BUILD_BINARY=ON BUILD_PYTHON=OFF BUILD_TEST=OFF ATEN_NO_TEST=ON python tools/build_libtorch.py --rerun-cmake
    _GLIBCXX_USE_CXX11_ABI=1 python tools/build_libtorch.py --rerun-cmake

    echo "Finished building libtorch in ${PWD}/build/"

fi

