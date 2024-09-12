
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

# If unset, assume it was installed into the current directory;
# see `install-build-deps-ubuntu.sh`.
CMAKE=${CMAKE:-"$(realpath ./cmake/bin/cmake)"}

git clone https://github.com/aws/aws-sdk-cpp.git --recurse-submodules
cd aws-sdk-cpp || exit
mkdir build
cd build || exit
cmake .. -DBUILD_ONLY="s3"
make
#sudo make install

