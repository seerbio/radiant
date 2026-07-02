#!/bin/bash

# "Strict mode"
set -eu -o pipefail

USE_BUILDX=${USE_BUILDX:-0}
DEBUG=${DEBUG:-1}
CACHE_FROM=${CACHE_FROM:-}
CACHE_TO=${CACHE_TO:-}
TEST_IMG=${TEST_IMG:-radiantdia-test:local}
APP_IMG=${APP_IMG:-radiantdia-app:local}

build_image() {
	local target=$1
	local image_ref=$2
	local cache_to=$3
	local -a docker_build_args

	if [[ "${USE_BUILDX}" == '1' ]]; then
		docker_build_args=(buildx build --load)
	else
		docker_build_args=(build)
	fi

	docker_build_args+=(--target "${target}" -t "${image_ref}")

	if [[ "${DEBUG}" == '0' ]]; then
		docker_build_args+=(-q)
	fi

	if [[ "${USE_BUILDX}" == '1' && -n "${CACHE_FROM}" ]]; then
		docker_build_args+=(--cache-from "${CACHE_FROM}")
	fi

	if [[ "${USE_BUILDX}" == '1' && -n "${cache_to}" ]]; then
		docker_build_args+=(--cache-to "${cache_to}")
	fi

	docker_build_args+=(.)
	docker "${docker_build_args[@]}"
}

echo "Building test container…"
build_image test "${TEST_IMG}" ""

# Run the default command from the `test` stage
echo "Running tests…"
docker run --rm "${TEST_IMG}"


# Build the final container. Note that this also builds
# the DEB package in a dedicated stage (but doesn't deploy
# it anywhere).
echo "Building app container…"
build_image app "${APP_IMG}" "${CACHE_TO}"

# Simply invoke the container with no arguments
# to ensure it's properly packaged
echo "Running container with no arguments…"
docker run --rm "${APP_IMG}"
