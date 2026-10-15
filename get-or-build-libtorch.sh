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

# Enable "strict mode" -- importantly this will cause the script to exit on any error
set -euo pipefail

# If unset, use `sudo apt-get` as the apt command
APT=${APT:-'sudo apt-get'}

# The version of PyTorch to use
PYTORCH_VERSION=${PYTORCH_VERSION:-'2.0.1'}

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
    wget "https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-${PYTORCH_VERSION}%2Bcpu.zip" -q -O "/tmp/libtorch-cxx11-abi-shared-with-deps-${PYTORCH_VERSION}%2Bcpu.zip"
    rm -rf /tmp/libtorch
    unzip "/tmp/libtorch-cxx11-abi-shared-with-deps-${PYTORCH_VERSION}%2Bcpu.zip" -d /tmp/
    rm "/tmp/libtorch-cxx11-abi-shared-with-deps-${PYTORCH_VERSION}%2Bcpu.zip"
    rm -rf "${PYTORCH_PREFIX_PATH}/pytorch/build/"
    mkdir -p "${PYTORCH_PREFIX_PATH}/pytorch/"
    mv /tmp/libtorch "${PYTORCH_PREFIX_PATH}/pytorch/build"

    echo "Done copying libtorch to ${PYTORCH_PREFIX_PATH}/pytorch/build/"

elif [[ "${ARCH}" == aarch64 || "${ARCH}" == arm64 ]] \
    && [[ "${PYTORCH_VERSION}" == 2.0.1 && "${PYTORCH_BUILD_FROM_SOURCE:-0}" != 1 ]]
then

    # The official ARM64 wheel contains the same C++ headers, CMake files and
    # libraries as LibTorch. Keep its adjacent torch.libs directory as well:
    # libtorch_cpu depends on the bundled ARM Compute and OpenMP libraries.
    echo "Using the pinned ARM64 LibTorch payload from the official PyTorch wheel..."
    ${APT} install --no-install-recommends -y python3 ca-certificates
    torch_wheel_path="/tmp/radiant-torch-2.0.1-aarch64.whl"
    wget -q \
        https://files.pythonhosted.org/packages/21/33/4925decd863ce88ed9190a4bd872b01c146243ee68db08c72923984fe335/torch-2.0.1-cp310-cp310-manylinux2014_aarch64.whl \
        -O "${torch_wheel_path}"
    printf '%s  %s\n' \
        359bfaad94d1cda02ab775dc1cc386d585712329bb47b8741607ef6ef4950747 \
        "${torch_wheel_path}" | sha256sum -c -
    rm -rf "${PYTORCH_PREFIX_PATH}/pytorch/build/"
    python3 - "${torch_wheel_path}" "${PYTORCH_PREFIX_PATH}/pytorch" <<'PY'
import shutil
import sys
import tempfile
from pathlib import Path
from zipfile import ZipFile

wheel, destination = map(Path, sys.argv[1:])
destination.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix="radiant-libtorch-") as temporary:
    extracted = Path(temporary)
    with ZipFile(wheel) as archive:
        archive.extractall(extracted)
    shutil.copytree(extracted / "torch", destination / "build")
    shutil.copytree(extracted / "torch.libs", destination / "torch.libs", dirs_exist_ok=True)
    for metadata in extracted.glob("*.dist-info"):
        shutil.copytree(metadata, destination / "build/share/licenses" / metadata.name,
                        dirs_exist_ok=True)
PY
    rm "${torch_wheel_path}"
    echo "Finished installing ARM64 LibTorch in ${PYTORCH_PREFIX_PATH}/pytorch/build/"

else

    echo "Found architecture: ${ARCH}"
    echo "Building libtorch..."

    # Build libtorch from sources
    mkdir -p "${PYTORCH_PREFIX_PATH}"
    cd "${PYTORCH_PREFIX_PATH}" || exit
    if [[ -d pytorch && ! -f pytorch/tools/build_libtorch.py ]]; then
        if [[ -f pytorch/build/share/cmake/Torch/TorchConfig.cmake ]]; then
            # Keep an earlier binary install when switching to source builds.
            # An extracted wheel is not a Git checkout, even though it occupies
            # the same installation prefix.
            torch_binary_backup=$(mktemp -d "${PWD}/pytorch-binary.XXXXXX")
            mv pytorch "${torch_binary_backup}/pytorch"
            echo "Preserved the previous binary installation in ${torch_binary_backup}/pytorch"
        else
            echo "Existing pytorch directory is neither a source checkout nor a recognized binary installation." >&2
            exit 2
        fi
    fi
    if [ ! -d pytorch ]; then git clone --depth 1 --recursive --shallow-submodules --branch "v${PYTORCH_VERSION}" https://github.com/pytorch/pytorch; fi
    cd pytorch || exit

    # Install Python, which is needed by the libtorch build
    ${APT} install -y python3.10 python-is-python3 python3-pip
    pip install setuptools pyyaml typing-extensions
    
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

