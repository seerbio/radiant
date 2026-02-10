#!/bin/bash

# "Strict mode"
set -eu -o pipefail

CONTAINER=${CONTAINER:-seer/radiant}

GIT_TAG_VERSION=$(git tag --points-at | grep -Po "^($CONTAINER/)?v?\K\d+\.\d+\.\d+.*$" | sort -Vr | head -n 1 || true)

export radiantdia_version=${radiantdia_version:-${GIT_TAG_VERSION}}

echo "Found radiantdia version '${radiantdia_version}'"

echo "Building deploy container…"
DEPLOY_IMG="${CONTAINER}-deb-build:${GIT_TAG_VERSION}"
docker build --target build-deb --build-arg radiantdia_version -t "${DEPLOY_IMG}" .

# Run the default command from the `deploy` stage.
# This builds the DEB and pushes it to S3.
echo "Running deployment…"
docker run --rm "${DEPLOY_IMG}"