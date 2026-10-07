#!/bin/bash

# "Strict mode"
set -eu -o pipefail

CONTAINER=${CONTAINER:-seer/radiant}
VERSION=${VERSION:-$(./get_release_version.sh)}
DEB_OUTPUT_DIR=${DEB_OUTPUT_DIR:?DEB_OUTPUT_DIR is required}
export radiantdia_version=${radiantdia_version:-${VERSION}}

echo "Found radiantdia version '${radiantdia_version}'"

echo "Building deploy container…"
DEPLOY_IMG="${CONTAINER}-deb-build:${radiantdia_version}"
docker_build_args=(
    build
    --load
    --target build-deb
    --build-arg "radiantdia_version=${radiantdia_version}"
    -t "${DEPLOY_IMG}"
)

docker buildx "${docker_build_args[@]}" .

mkdir -p "${DEB_OUTPUT_DIR}"
HOST_DEB_OUTPUT_DIR=$(realpath "${DEB_OUTPUT_DIR}")
HOST_UID=$(id -u)
HOST_GID=$(id -g)

echo "Exporting DEB artifacts…"
docker run --rm \
    --user "${HOST_UID}:${HOST_GID}" \
    -v "${HOST_DEB_OUTPUT_DIR}:/tmp/deb-out" \
    -e "LOCAL_DEB_OUTPUT_DIR=/tmp/deb-out" \
    "${DEPLOY_IMG}"
