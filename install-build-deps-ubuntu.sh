#!/usr/bin/bash

# install-build-deps-ubuntu.sh
# 
# Set up dependencies for building from scratch on Ubuntu

# If unset, use `sudo apt-get` as the apt command
APT=${APT:-'sudo apt-get'}

# If unset, install `cmake` into the current directory
CMAKE_PREFIX=${CMAKE_PREFIX:-'./cmake'}

# Get initial requirements
${APT} install --no-install-recommends -y ca-certificates lsb-release wget

wget https://apache.jfrog.io/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
${APT} install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
${APT} update

${APT} install --no-install-recommends -y \
    git \
    build-essential \
    unzip
    # Possibly not required; artifacts from libtorch build development?
    #python3.10 \
    #python-is-python3 \
    #python3-pip \
    libboost-all-dev \
    qtbase5-dev \
    libarrow-dev \
    libparquet-dev \

# Possibly not required; artifact from libtorch build development?
#pip install --no-cache-dir setuptools pyyaml

# Install latest CMAKE > 3.17
wget https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2-Linux-`uname -m`.sh -q -O /tmp/cmake-install.sh
chmod u+x /tmp/cmake-install.sh
mkdir /usr/bin/cmake
/tmp/cmake-install.sh --skip-license --prefix=${CMAKE_PREFIX}
rm /tmp/cmake-install.sh
