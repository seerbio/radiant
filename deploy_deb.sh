#!/bin/bash

# "Strict mode"
set -eu -o pipefail

CONTAINER=${CONTAINER:-seer/radiant}
VERSION=${VERSION:-$(./get_release_version.sh)}
export radiantdia_version=${radiantdia_version:-${VERSION}}

echo "Found radiantdia version '${radiantdia_version}'"

echo "Building deploy container…"
DEPLOY_IMG="${CONTAINER}-deb-build:${radiantdia_version}"
docker build --target build-deb --build-arg radiantdia_version -t "${DEPLOY_IMG}" .

# Run the default command from the `deploy` stage.
# This builds the DEB and pushes it to S3.
echo "Running deployment…"
docker_run_args=(--rm)

if [[ -n "${DEB_OUTPUT_DIR:-}" ]]; then
    mkdir -p "${DEB_OUTPUT_DIR}"
    HOST_DEB_OUTPUT_DIR=$(realpath "${DEB_OUTPUT_DIR}")
    HOST_UID=$(id -u)
    HOST_GID=$(id -g)

    echo "Will also copy the DEB to '${HOST_DEB_OUTPUT_DIR}'"
    echo "Will run the uploader as uid:gid '${HOST_UID}:${HOST_GID}'"

    docker_run_args+=(
        --user "${HOST_UID}:${HOST_GID}"
        -v "${HOST_DEB_OUTPUT_DIR}:/tmp/deb-out"
        -e "HOME=/tmp"
        -e "LOCAL_DEB_OUTPUT_DIR=/tmp/deb-out"
    )
fi

docker run "${docker_run_args[@]}" "${DEPLOY_IMG}"
