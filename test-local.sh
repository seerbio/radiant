#!/usr/bin/bash

# build-local.sh
#
# Build Pythia in the current directory.
#
# FOR DEVELOPMENT PURPOSES ONLY!!
#
# Assumes that dependencies were installed using the
# default configuration of the install/get/build scripts.
# For more information see `CONTRIBUTING.md`.

# "Strict mode"
set -euo pipefail

cmake/bin/ctest --test-dir build/ --output-on-failure
