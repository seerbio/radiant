#!/usr/bin/bash

# install-build-deps-ubuntu.sh
# 
# Set up dependencies for building from scratch on Ubuntu

# If unset, use `sudo apt-get` as the apt command
APT=${APT:-'sudo apt-get'}

# The version of CMake to install; must be > 3.17!
CMAKE_VERSION=${CMAKE_VERSION:-'3.23.2'}

# If unset, install `cmake` into the current directory
CMAKE_PREFIX=${CMAKE_PREFIX:-'./cmake'}

# Enable "strict mode" -- importantly this will cause the script to exit on any error
set -euo pipefail

# Get initial requirements
${APT} install --no-install-recommends -y ca-certificates lsb-release wget

# Allow overriding the distribution name -- on Pop!OS run with `DISTRO=`.
# If DISTRO is unset, use `lsb_release` to determine the distribution name.
# If DISTRO is set to an empty value or `lsb_release` returns nothing, fall back to the default 'ubuntu'.
DISTRO=${DISTRO-"$(lsb_release --id --short | tr '[:upper:]' '[:lower:]')"}
DISTRO=${DISTRO:-ubuntu}

CODENAME=$(lsb_release --codename --short)

wget -P /tmp/ "https://apache.jfrog.io/artifactory/arrow/${DISTRO}/apache-arrow-apt-source-latest-${CODENAME}.deb"
${APT} install -y -V "/tmp/apache-arrow-apt-source-latest-${CODENAME}.deb"
rm "/tmp/apache-arrow-apt-source-latest-${CODENAME}.deb"
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
    libxml2-dev \
    libqt5opengl5-dev \
    libqt5svg5-dev

# Install latest CMAKE > 3.17
wget "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-$(uname -m).sh" -q -O /tmp/cmake-install.sh
chmod u+x /tmp/cmake-install.sh
mkdir -p "${CMAKE_PREFIX}"
/tmp/cmake-install.sh --skip-license --prefix="${CMAKE_PREFIX}"
rm /tmp/cmake-install.sh
