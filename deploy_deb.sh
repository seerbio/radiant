#!/bin/bash

# "Strict mode"
set -eu -o pipefail

CONTAINER=${CONTAINER:-seer/pythia-dia}

GIT_TAG_VERSION=$(git tag --points-at | grep -Po "^($CONTAINER/)?v?\K\d+\.\d+\.\d+.*$" | sort -Vr | head -n 1 || true)

export pythiadia_version=${pythiadia_version:-${GIT_TAG_VERSION}}

echo "Found pythiadia version '${pythiadia_version}'"

echo "Building deploy container…"
DEPLOY_IMG="${CONTAINER}-deb-build:${GIT_TAG_VERSION}"
docker build --target build-deb --build-arg pythiadia_version -t "${DEPLOY_IMG}" .

# Run the default command from the `deploy` stage.
# This builds the DEB and pushes it to S3.
echo "Running deployment…"
docker run --rm "${DEPLOY_IMG}"