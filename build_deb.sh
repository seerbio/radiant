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

radiantdia_version="${radiantdia_version:-0.0-dev}"

package_dir="${package_dir}+${ARCH}"

radiant_bin="${radiant_bin:=${package_dir}/usr/local/bin/radiant}"
radiant_lib="${radiant_lib:=${package_dir}/usr/lib}"

mkdir -p "${radiant_bin}"
mkdir -p "${radiant_lib}"

echo "Copying DEB contents to package_dir…"

#cp bin/IstrosLibraryBuilder "$radiant_bin/IstrosLibraryBuilder"
#cp bin/LycaonMzMLTransmorgitron "$radiant_bin/LycaonMzMLTransmorgitron"
#cp bin/OrpheusFastaShredder "$radiant_bin/OrpheusFastaShredder"
cp bin/RadiantDIA "$radiant_bin/RadiantDIA"
cp bin/FeatureOptimizerPrime "$radiant_bin/FeatureOptimizerPrime"
cp AlgorithmsFFLib/libAlgorithmsFFLib.so "$radiant_lib/libAlgorithmsFFLib.so"
cp ChemLib/libChemLib.so "$radiant_lib/libChemLib.so"
cp EigenLib/libEigenLib.so "$radiant_lib/libEigenLib.so"
cp FileReadersLib/libFileReadersLib.so "$radiant_lib/libFileReadersLib.so"
cp UtilsLib/libUtilsLib.so "$radiant_lib/libUtilsLib.so"
cp WorkFlowsLib/libWorkFlowsLib.so "$radiant_lib/libWorkFlowsLib.so"
cp PyTorchLib/libPyTorchLib.so "$radiant_lib/libPyTorchLib.so"

cp /src/pytorch/build/lib/libtorch.so "$radiant_lib/"
cp /src/pytorch/build/lib/libtorch_cpu.so "$radiant_lib/"
cp /src/pytorch/build/lib/libc10.so "$radiant_lib/"

# Prebuilt libtorch bundles a libgomp runtime; source builds do not.
for f in /src/pytorch/build/lib/libgomp-*.so*; do
    [ -e "$f" ] || break
    cp "$f" "$radiant_lib/"
done

# Create symlinks
# Must match ${radiant_bin} above, but relative to the target!
ln -s ../../../usr/local/bin/radiant/RadiantDIA "${package_dir}/usr/local/bin/RadiantDIA"

mkdir -p "${package_dir}/DEBIAN/"
cp "control.${ARCH}" "${package_dir}/DEBIAN/control"
echo "Version: ${radiantdia_version}" >> "${package_dir}/DEBIAN/control"

echo "Building RadiantDIA DEB…"

dpkg-deb --build "$package_dir"
