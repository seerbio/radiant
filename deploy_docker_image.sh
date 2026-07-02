#!/bin/bash

# "Strict mode"
set -eu -o pipefail

CONTAINER=${CONTAINER:-seer/radiant}
REPOSITORY=${REPOSITORY:-docker.io}
TAG_SUFFIX=${TAG_SUFFIX:?TAG_SUFFIX is required}
CACHE_SCOPE=${CACHE_SCOPE:-}

VERSION=${VERSION:-$(./get_release_version.sh)}

echo "Found RadiantDIA version '${VERSION}'"

TAGS=("${VERSION}${TAG_SUFFIX}")
if [[ -n "${1:-}" ]]; then
  read -r -a EXTRA_TAGS <<< "${1}"
  TAGS+=("${EXTRA_TAGS[@]}")
fi

docker_build_args=(
  build
  --load
  --target app
  --build-arg "radiantdia_version=${VERSION}"
)

if [[ -n "${CACHE_SCOPE}" ]]; then
  docker_build_args+=(
    --cache-from "type=gha,scope=${CACHE_SCOPE}"
    --cache-to "type=gha,scope=${CACHE_SCOPE},mode=max"
  )
fi

for TAG in "${TAGS[@]}"; do
  docker_build_args+=(-t "${REPOSITORY}/${CONTAINER}:${TAG}")
done

echo "Building container image(s): ${TAGS[*]}"
docker buildx "${docker_build_args[@]}" .

for TAG in "${TAGS[@]}"; do
  echo "Pushing ${REPOSITORY}/${CONTAINER}:${TAG}"
  docker push "${REPOSITORY}/${CONTAINER}:${TAG}"
done
