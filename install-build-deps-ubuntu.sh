#!/usr/bin/bash

# install-build-deps-ubuntu.sh
# 
# Set up dependencies for building from scratch on Ubuntu

# If unset, use `sudo apt-get` as the apt command
APT=${APT:-'sudo apt-get'}

# If unset, install `cmake` into the current directory
CMAKE_PREFIX=${CMAKE_PREFIX:-'./cmake'}

# Enable "strict mode" -- importantly this will cause the script to exit on any error
set -euo pipefail

# Get initial requirements
${APT} install --no-install-recommends -y ca-certificates lsb-release wget

# Allow overriding the distribution name -- on Pop!OS run with `DISTRO=`
# Fallback to using lsb_release only occurs when DISTRO is not set; if blank we will fall back to default below.
DISTRO=${DISTRO-"$(lsb_release --id --short | tr '[:upper:]' '[:lower:]')"}
CODENAME=$(lsb_release --codename --short)

wget -P /tmp/ https://apache.jfrog.io/artifactory/arrow/${DISTRO:-ubuntu}/apache-arrow-apt-source-latest-${CODENAME}.deb
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
wget "https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2-Linux-$(uname -m).sh" -q -O /tmp/cmake-install.sh
chmod u+x /tmp/cmake-install.sh
mkdir -p "${CMAKE_PREFIX}"
/tmp/cmake-install.sh --skip-license --prefix="${CMAKE_PREFIX}"
rm /tmp/cmake-install.sh
