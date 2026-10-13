#!/bin/bash

# "Strict mode"
set -eu -o pipefail

CONTAINER=${CONTAINER:-seer/radiant}
export VERSION=${VERSION:-$(./get_release_version.sh)}

echo "Found RadiantDIA version '${VERSION}'"

MAJVER=$(grep -Po '^.+?(?=\..+)' <<< "$VERSION") \
  || (echo "Could not parse major version from $VERSION" && exit 1)
MINVER=$(grep -Po '^.+?\..+?(?=\..+)' <<< "$VERSION") \
  || (echo "Could not parse minor version from $VERSION" && exit 1)

# Build and push a manifest with the semantic version.
echo "Building multi-arch manifest ${REPOSITORY}/${CONTAINER}:${VERSION}"
docker manifest create \
  "${REPOSITORY}/${CONTAINER}:${VERSION}" \
  "${REPOSITORY}/${CONTAINER}:${VERSION}-x86_64" \
  "${REPOSITORY}/${CONTAINER}:${VERSION}-aarch64"
docker manifest push --purge "${REPOSITORY}/${CONTAINER}:${VERSION}"

# Build/push manifest for minor/major versions
for TAG in \
  "${MINVER}" \
  "${MAJVER}" ;
do
  echo "Tagging ${REPOSITORY}/${CONTAINER}:${VERSION} as ${CONTAINER}:${TAG}"

  docker manifest create \
    "${REPOSITORY}/${CONTAINER}:${TAG}" \
    "${REPOSITORY}/${CONTAINER}:${VERSION}-x86_64" \
    "${REPOSITORY}/${CONTAINER}:${VERSION}-aarch64"
  docker manifest push --purge "${REPOSITORY}/${CONTAINER}:${TAG}"
done

# If a second argument was given, treat it as a single string of space-separated additional tags.
# Build/push manifest with these tag(s)
for TAG in $(echo "${1:-}") ; do
  echo "Tagging ${REPOSITORY}/${CONTAINER}:${VERSION} as ${CONTAINER}:${TAG}"

  docker manifest create \
    "${REPOSITORY}/${CONTAINER}:${TAG}" \
    "${REPOSITORY}/${CONTAINER}:${VERSION}-x86_64" \
    "${REPOSITORY}/${CONTAINER}:${VERSION}-aarch64"
  docker manifest push --purge "${REPOSITORY}/${CONTAINER}:${TAG}"
done
