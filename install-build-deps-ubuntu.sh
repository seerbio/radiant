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

wget -P /tmp/ "https://apache.jfrog.io/artifactory/arrow/$(lsb_release --id --short | tr '[:upper:]' '[:lower:]')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb"
${APT} install -y -V "/tmp/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb"
rm "/tmp/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb"
${APT} update

${APT} install --no-install-recommends -y \
    git \
    build-essential \
    unzip \
    libboost-all-dev \
    qtbase5-dev \
    libarrow-dev \
    libparquet-dev \
    libseqan2-dev \
    libsqlite3-dev \
    libxml2-dev

# Install latest CMAKE > 3.17
wget "https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2-Linux-$(uname -m).sh" -q -O /tmp/cmake-install.sh
chmod u+x /tmp/cmake-install.sh
mkdir -p "${CMAKE_PREFIX}"
/tmp/cmake-install.sh --skip-license --prefix="${CMAKE_PREFIX}"
rm /tmp/cmake-install.sh
