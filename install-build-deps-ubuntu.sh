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

wget -P /tmp/ "https://apache.jfrog.io/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb"
${APT} install -y -V "/tmp/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb"
rm "/tmp/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb"
${APT} update

${APT} install --no-install-recommends -y \
    build-essential \
    qtbase5-dev \
    libseqan2-dev \
    git \
    unzip \
    zip \
    curl \
    tar \
    pkg-config \
    flex \
    bison \
    autoconf \
    automake \
    libtool \
    python-is-python3 \
    m4 \
    gettext \
    libboost-all-dev \
    qtbase5-dev \
    python3.10  \
    python3.10-venv \
    python3-pip  \
    pkg-config  \
    libcurl4-openssl-dev \
    libssl-dev \
    unzip \
    zip \
    git \
    flex \
    uuid-dev \
    zlib1g-dev \
    autoconf \
    automake \
    autoconf-archive \
    libtool \
    bison \
    libxi-dev \
    libxtst-dev \
    libx11-dev \
    libxft-dev \
    libxext-dev \
    pkg-config \
    ninja-build \
    libglib2.0-dev \
    libgdk-pixbuf2.0-dev \
    libpango1.0-dev \
    libcairo2-dev \
    libxrandr-dev \
    gfortran \
    libgirepository1.0-dev \
    m4 \
    gettext \
    libgoogle-glog-dev \
    libarchive-dev \
    libxkbcommon-dev \
    libpcre2-dev \
    libglib2.0-dev \
    libdbus-1-dev \
    python3-jinja2 \
    libgflags-dev

# Install latest CMAKE > 3.23
wget "https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2-Linux-$(uname -m).sh" -q -O /tmp/cmake-install.sh
chmod u+x /tmp/cmake-install.sh
mkdir -p "${CMAKE_PREFIX}"
/tmp/cmake-install.sh --skip-license --prefix="${CMAKE_PREFIX}"
rm /tmp/cmake-install.sh

git clone https://github.com/microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg update
./vcpkg/vcpkg install \
    boost \
    arrow \
    parquet \
    curl \
    openssl \
    thrift \
    gflags \
    glog \
    libarchive \
    pcre2 \
    dbus \
    aws-sdk-cpp \
    zstd


./vcpkg/vcpkg integrate install

#./vcpkg/vcpkg install aws-sdk-cpp
#./vcpkg/vcpkg install utf8proc
##./vcpkg/vcpkg install libtorch
##./vcpkg/vcpkg install lz4
##./vcpkg/vcpkg install zstd
##./vcpkg/vcpkg install re2
##./vcpkg/vcpkg install parquet
##./vcpkg/vcpkg install arrow
#./vcpkg/vcpkg integrate install
