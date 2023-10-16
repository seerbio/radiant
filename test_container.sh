#!/bin/bash

# "Strict mode"
set -eu -o pipefail

echo "Building test container…"

DEBUG=${DEBUG:-1}

if [ ${DEBUG} != '0' ];
then
	# Build the container, with log output
	docker build --target test .
fi

# Re-run the build command to get the ID (should be cached)
TEST_IMG=$(docker build -q --target test .)

# Run the default command from the `test` stage
echo "Running tests…"
docker run --rm "${TEST_IMG}"


# Build the final container. Note that this also builds
# the DEB package in a dedicated stage (but doesn't deploy
# it anywhere).
echo "Building app container…"

if [ ${DEBUG} != '0' ];
then
	# Build the container, with log output
	docker build --target app .
fi

IMG=$(docker build -q --target app .)

# Simply invoke the container with no arguments
# to ensure it's properly packaged
echo "Running container with no arguments…"
docker run --rm "${IMG}"
