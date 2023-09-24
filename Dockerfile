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
    && apt-get upgrade -y \
    && apt-get install --no-install-recommends -y ca-certificates wget build-essential  qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools hdf5-tools  \
    libcurl4-openssl-dev  libhdf5-dev libbrotli-dev libboost-all-dev libutf8proc2 libre2-5 libsnappy1v5 libthrift-0.13.0 \
    unzip=6.0-25ubuntu1.1 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install latest CMAKE > 3.17
RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2-Linux-x86_64.sh -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /usr/bin/cmake \
      && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
      && rm /tmp/cmake-install.sh

ENV PATH="/usr/bin/cmake/bin:${PATH}"

# Copy project source into the container
COPY ./ /src/PythiaDIACpp/

# https://pytorch.org
RUN wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip -q -O ./libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
RUN unzip ./libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip

# Build the project in /app/
#WORKDIR /app/
#RUN cmake -DCMAKE_BUILD_TYPE=Release -S /src/PythiaDIACpp/ -B /app/ \
#    && make
#
#################################################
##
## Test stage
##
## Here we add test dependencies to the build container
## and define an ENTRYPOINT that will run them.
##
#################################################
#FROM build AS test
#
##
## This entrypoint allows running tests (with coverage)
## by building with `--target test` and then running the
## resulting container. See README.md for more.
##
#WORKDIR /app/
#CMD ["ctest"]
#
#
################################################
##
## Deploy stage
##
## Here we put everything in its right place for
## Debian Package Deployment and build the DEB.
## Once built, deployment is as easy as running
## this stage's default command:
##
##     $ docker run --rm -it $(docker build --target deploy .)
##
#################################################
#FROM build AS deploy
#
## Install Python and dependencies
#RUN apt-get update \
#    && apt-get upgrade -y \
#    && apt-get install --no-install-recommends -y python3.9 python-is-python3 python3-pip \
#    && apt-get autoremove -y \
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/* \
#    && pip install --no-cache-dir awscli boto3
#
## Copy some extra contents into /app/
#RUN cp \
#    /src/PythiaDIACpp/control \
#    /src/PythiaDIACpp/build_deb.sh \
#    /src/PythiaDIACpp/s3_package_uploader.py \
#    /app/
#
#RUN cp /src/PythiaDIACpp/ThirdPartyLibs/arrow_parquet/release/libarrow.so.1100 /usr/lib/libarrow.so.1100
#RUN cp /src/PythiaDIACpp/ThirdPartyLibs/arrow_parquet/release/libparquet.so.1100 /usr/lib/libparquet.so.1100
#
#WORKDIR /app/
#
#ARG pythia_dia_version=0.1
#ENV package_dir=pythia_dia_${pythia_dia_version}
#ENV PACKAGE_NAME=${package_dir}.deb
#
## This should work with the default entrypoint to build and deploy.
#CMD ["/bin/bash -c ./build_deb.sh && python s3_package_uploader.py"]
#
##################################################
###
### App stage
###
### Here we build the final container used to run the app in production.
### Must be the last stage for compatibility with GitHub Actions build.
###
##################################################
##FROM base AS app
##
###
### Set labels
###
##LABEL author="Seer, Inc."
##LABEL description="PythiaDIACpp"
##
##RUN apt-get update \
##    && apt-get install --no-install-recommends -y build-essential cmake qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools  \
##    hdf5-tools libhdf5-dev libboost-all-dev libutf8proc2 libre2-5 libsnappy1v5 libthrift-0.13.0 libcurl4-openssl-dev  libbrotli-dev \
##    && apt-get autoremove -y \
##    && apt-get clean \
##    && rm -rf /var/lib/apt/lists/*
##
##COPY --from=deploy /app/ /app/
###COPY --from=build /app/AlgorithmsLib/ /app/AlgorithmsLib/
###COPY --from=build /app/ChemLib/ /app/ChemLib/
###COPY --from=build /app/EigenLib/ /app/EigenLib/
###COPY --from=build /app/FileReadersLib/ /app/FileReadersLib/
###COPY --from=build /app/MachineLrnLib/ /app/MachineLrnLib/
###COPY --from=build /app/UtilsLib/ /app/UtilsLib/
###COPY --from=build /app/WorkFlowsLib/ /app/WorkFlowsLib/
###COPY --from=build /app/ThirdPartyLibs/arrow_parquet/release/. /app/bin/
##
### Set up a dedicated folder as the normal working dir
##WORKDIR /work/
##
### Using this entrypoint means the "command" passed to `docker run` will be arguments to
### this binary (e.g. `docker run seer/pythia -h`). To run a different binary requires
### overriding the entrypoint (e.g. `docker run --entrypoint /app/bin/Pythia)
##ENTRYPOINT ["/app/bin/PythiaDIACpp"]
