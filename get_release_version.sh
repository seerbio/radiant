#!/bin/bash

# "Strict mode"
set -eu -o pipefail

TAG_MATCH_CONTAINER=${TAG_MATCH_CONTAINER:-${CONTAINER:-}}

if [[ -n "${TAG_MATCH_CONTAINER}" ]]; then
  TAG_PREFIX_REGEX="(${TAG_MATCH_CONTAINER}/)?"
else
  TAG_PREFIX_REGEX=""
fi

VERSION=$(git tag --points-at | grep -Po "^${TAG_PREFIX_REGEX}v?\K\d+\.\d+\.\d+.*$" | sort -Vr | head -n 1 || true)

if [[ -z "${VERSION}" ]]; then
  if [[ -n "${TAG_MATCH_CONTAINER}" ]]; then
    echo "Could not determine VERSION from git tags matching '${TAG_MATCH_CONTAINER}'" >&2
  else
    echo "Could not determine VERSION from git tags on HEAD" >&2
  fi
  exit 1
fi

printf '%s\n' "${VERSION}"
