#!/bin/bash

# "Strict mode"
set -eu -o pipefail

export pythia_version=${pythia_version}

echo "Building deploy container…"
DEPLOY_IMG=$(docker build -q --target deploy --build-arg pythia_dia_version .)

# Run the default command from the `deploy` stage.
# This builds the DEB and pushes it to S3.
echo "Running deployment…"
docker run --rm "${DEPLOY_IMG}"