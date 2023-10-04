###############################################
#
# Base stage
#
# Here we set up common configuration and runtime dependencies
# that are shared between the build stage and final image.
# At this stage we also copy the relevant app files so they're
# accessible in all later stages.
#
################################################

# Based on Ubuntu 20.04 LTS ("Focal Fossa")
FROM ubuntu:20.04 AS base

#
# Set locales to UTF-8
#
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

#
# Set timezone
#
ENV TZ=US/Pacific
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

################################################
#
# Build stage
#
# Here we install development/build tools and fetch dependencies.
#
################################################


# ###
# # DREW'S MANUAL BUILD STEPS:
# ###
#
# https://hub.docker.com/r/arrowdev/amd64-debian-10-cpp
# docker pull arrowdev/amd64-debian-10-cpp
# sudo docker run -it --entrypoint /bin/bash arrowdev/amd64-debian-10-cpp
# git clone https://github.com/apache/arrow.git
# cd arrow/cpp
# mkdir build
# cd build
# cmake .. -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DARROW_WITH_ZLIB=ON -DPARQUET_BUILD_STATIC=ON -DARROW_PARQUET=ON -DARROW_WITH_SNAPPY=ON -DPARQUET_BUILD_EXECUTABLES=ON
# make

### INSTRUCTIONS TO PULL FILES FROM DOCKER CONTAINER
# zip the release directory

# Open separate termial window
# sudo docker ps
# sudo docker exec <container_id> chmod 777 /path/to/zipped/directory

# docker cp <container_name>:path/to/file -
# docker cp <container_name>:path/to/file - > file.txt

# sudo docker cp 7449c6433401:/arrow/cpp/arrow.zip -
# sudo docker cp 7449c6433401:/arrow/cpp/arrow.zip - > arrow.zip

# sudo docker cp d99e58c57b48:/arrow/cpp/build/arrow_parquet.zip -
# sudo docker cp d99e58c57b48:/arrow/cpp/build/arrow_parquet.zip - > arrow_parquet.zip
########################################################


FROM apache/arrow-dev:arm64v8-ubuntu-20.04-cpp AS build-arrow

WORKDIR /src/
RUN git clone https://github.com/apache/arrow.git \
    && mkdir arrow/cpp/build/ \
    && cmake -S arrow/cpp/ -B arrow/cpp/build/ -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DARROW_WITH_ZLIB=ON -DPARQUET_BUILD_STATIC=ON -DARROW_PARQUET=ON -DARROW_WITH_SNAPPY=ON -DPARQUET_BUILD_EXECUTABLES=ON \
    && cd arrow/cpp/build \
    && make -j

FROM base AS build

#
# Update, upgrade packages. Install build dependencies.
#
# We use a single RUN statement to reduce the number of image layers
# created by docker, which reduces overall image size.
#
# After package installation, we remove any unnecessary packages
# and run `apt-get clean` to remove any downloaded package archives.
#
RUN apt-get update \
#    && apt-get install -y ca-certificates lsb-release wget \
#    && wget https://apache.jfrog.io/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb \
#    && apt-get install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb \
    && apt-get update \
    && apt-get upgrade -y \
    && apt-get install --no-install-recommends -y ca-certificates wget build-essential  qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools hdf5-tools  \
        libcurl4-openssl-dev  libhdf5-dev libbrotli-dev libboost-all-dev libutf8proc2 libre2-5 libsnappy1v5 libthrift-0.13.0 \
#        libarrow-dev libparquet-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install latest CMAKE > 3.17
RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2-Linux-`uname -m`.sh -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /usr/bin/cmake \
      && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
      && rm /tmp/cmake-install.sh

ENV PATH="/usr/bin/cmake/bin:${PATH}"

# Copy Arrow/Parquet libraries
COPY --from=build-arrow /src/arrow/cpp/build/release /src/arrow/cpp/build/release

# Copy project source into the container
COPY ./ /src/PythiaDIACpp/

# Build the project in /app/
WORKDIR /app/
RUN cmake -S /src/PythiaDIACpp/ -B /app/ -DCMAKE_PREFIX_PATH=/src/arrow/cpp/build/ -DCMAKE_BUILD_TYPE=Release \
    && make -j

################################################
#
# Test stage
#
# Here we add test dependencies to the build container
# and define an ENTRYPOINT that will run them.
#
################################################
FROM build AS test

#
# This entrypoint allows running tests (with coverage)
# by building with `--target test` and then running the
# resulting container. See README.md for more.
#
WORKDIR /app/
CMD ["ctest"]


###############################################
#
# Deploy stage
#
# Here we put everything in its right place for
# Debian Package Deployment and build the DEB.
# Once built, deployment is as easy as running
# this stage's default command:
#
#     $ docker run --rm -it $(docker build --target deploy .)
#
################################################
FROM build AS deploy

# Install Python and dependencies
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install --no-install-recommends -y python3.9 python-is-python3 python3-pip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir awscli boto3

# Copy some extra contents into /app/
RUN cp \
    /src/PythiaDIACpp/control \
    /src/PythiaDIACpp/build_deb.sh \
    /src/PythiaDIACpp/s3_package_uploader.py \
    /app/

RUN cp /src/PythiaDIACpp/ThirdPartyLibs/arrow_parquet/release/libarrow.so.1100 /usr/lib/libarrow.so.1100
RUN cp /src/PythiaDIACpp/ThirdPartyLibs/arrow_parquet/release/libparquet.so.1100 /usr/lib/libparquet.so.1100

WORKDIR /app/

ARG pythia_dia_version=0.1
ENV package_dir=pythia_dia_${pythia_dia_version}
ENV PACKAGE_NAME=${package_dir}.deb

# This should work with the default entrypoint to build and deploy.
CMD ["/bin/bash -c ./build_deb.sh && python s3_package_uploader.py"]

#################################################
##
## App stage
##
## Here we build the final container used to run the app in production.
## Must be the last stage for compatibility with GitHub Actions build.
##
#################################################
#FROM base AS app
#
##
## Set labels
##
#LABEL author="Seer, Inc."
#LABEL description="PythiaDIACpp"
#
#RUN apt-get update \
#    && apt-get install --no-install-recommends -y build-essential cmake qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools  \
#    hdf5-tools libhdf5-dev libboost-all-dev libutf8proc2 libre2-5 libsnappy1v5 libthrift-0.13.0 libcurl4-openssl-dev  libbrotli-dev \
#    && apt-get autoremove -y \
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/*
#
#COPY --from=deploy /app/ /app/
##COPY --from=build /app/AlgorithmsLib/ /app/AlgorithmsLib/
##COPY --from=build /app/ChemLib/ /app/ChemLib/
##COPY --from=build /app/EigenLib/ /app/EigenLib/
##COPY --from=build /app/FileReadersLib/ /app/FileReadersLib/
##COPY --from=build /app/MachineLrnLib/ /app/MachineLrnLib/
##COPY --from=build /app/UtilsLib/ /app/UtilsLib/
##COPY --from=build /app/WorkFlowsLib/ /app/WorkFlowsLib/
##COPY --from=build /app/ThirdPartyLibs/arrow_parquet/release/. /app/bin/
#
## Set up a dedicated folder as the normal working dir
#WORKDIR /work/
#
## Using this entrypoint means the "command" passed to `docker run` will be arguments to
## this binary (e.g. `docker run seer/pythia -h`). To run a different binary requires
## overriding the entrypoint (e.g. `docker run --entrypoint /app/bin/Pythia)
#ENTRYPOINT ["/app/bin/PythiaDIACpp"]
