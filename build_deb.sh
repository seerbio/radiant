#!/bin/bash
#
# build_deb.sh
#
# Run from a directory containing an already-built version of Pythia,
# e.g. the `/app/` directory of the container built by `Dockerfile`.
#
# You MUST define `$package_dir`, which is the location where the
# contents of the DEB will be assembled. They will not be cleaned up.
#
# The created DEB should be located at `${package_dir}-${ARCH}.deb`,
# where `${ARCH}` is (typically) `amd64` or `arm64`.

# "Strict mode"
set -eu -o pipefail

# Architecture value used for the DEB, typically `amd64` or `arm64`
# (Note: this differs from uname-based strings used elsewhere).
# This string will be appended to the package name.
ARCH=$(dpkg-architecture | grep 'DEB_BUILD_ARCH=' | cut -d = -f 2)

# Required argument (as env var)
package_dir="${package_dir}+${ARCH}"

pythia_bin="${pythia_bin:=${package_dir}/usr/local/bin/PythiaDIACpp}"
pythia_lib="${pythia_lib:=${package_dir}/usr/lib}"

mkdir -p "${pythia_bin}"
mkdir -p "${pythia_lib}"

echo "Copying DEB contents to package_dir…"

#cp bin/IstrosLibraryBuilder "$pythia_bin/IstrosLibraryBuilder"
#cp bin/LycaonMzMLTransmorgitron "$pythia_bin/LycaonMzMLTransmorgitron"
#cp bin/OrpheusFastaShredder "$pythia_bin/OrpheusFastaShredder"
cp bin/PythiaDIA "$pythia_bin/PythiaDIA"
cp bin/iRT_Model.json "$pythia_bin/iRT_Model.json"
cp bin/MS2_Mono_Model.json "$pythia_bin/MS2_Mono_Model.json"
cp bin/MS2_Charge_Model.json "$pythia_bin/MS2_Charge_Model.json"
cp AlgorithmsFFLib/libAlgorithmsFFLib.so "$pythia_lib/libAlgorithmsFFLib.so"
cp ChemLib/libChemLib.so "$pythia_lib/libChemLib.so"
cp EigenLib/libEigenLib.so "$pythia_lib/libEigenLib.so"
cp FileReadersLib/libFileReadersLib.so "$pythia_lib/libFileReadersLib.so"
cp UtilsLib/libUtilsLib.so "$pythia_lib/libUtilsLib.so"
cp WorkFlowsLib/libWorkFlowsLib.so "$pythia_lib/libWorkFlowsLib.so"
cp PyTorchLib/libPyTorchLib.so "$pythia_lib/libPyTorchLib.so"

cp /src/pytorch/build/lib/* "$pythia_lib/"

mkdir -p "${package_dir}/DEBIAN/"
cp "control.${ARCH}" "${package_dir}/DEBIAN/control"

echo "Building PythiaDIA DEB…"

dpkg-deb --build "$package_dir"
