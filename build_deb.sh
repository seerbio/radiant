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
# The created DEB is located at `${package_dir}+${ARCH}.deb`,
# where `${ARCH}` is (typically) `amd64` or `arm64`.
# TORCH_LIBRARY_DIR can select an unpacked libtorch or Python wheel's torch/lib.
# THERMO_RAW_READER_DIR can select a self-contained publish outside bin/thermoraw.

# "Strict mode"
set -eu -o pipefail

# Architecture value used for the DEB, typically `amd64` or `arm64`
# (Note: this differs from uname-based strings used elsewhere).
# This string will be appended to the package name.
ARCH=$(dpkg-architecture -qDEB_HOST_ARCH)
case "$ARCH" in
    amd64) elf_machine='Advanced Micro Devices X86-64' ;;
    arm64) elf_machine='AArch64' ;;
    *) printf 'Unsupported package architecture: %s\n' "$ARCH" >&2; exit 2 ;;
esac

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

# Keep the complete published helper beside RadiantDIA for automatic discovery.
# Shipping only its apphost would silently require a separately installed .NET.
THERMO_RAW_READER_DIR="${THERMO_RAW_READER_DIR:-bin/thermoraw}"
for file in RadiantThermoRawReader RadiantThermoRawReader.dll \
    RadiantThermoRawReader.runtimeconfig.json libhostpolicy.so libcoreclr.so; do
    if [ ! -f "${THERMO_RAW_READER_DIR}/${file}" ]; then
        printf 'Missing self-contained Thermo RAW reader file: %s/%s\n' \
            "$THERMO_RAW_READER_DIR" "$file" >&2
        exit 2
    fi
done
if [ ! -x "${THERMO_RAW_READER_DIR}/RadiantThermoRawReader" ] \
    || ! grep -q '"includedFrameworks"' "${THERMO_RAW_READER_DIR}/RadiantThermoRawReader.runtimeconfig.json"; then
    printf 'Thermo RAW reader must be an executable, self-contained .NET publish.\n' >&2
    exit 2
fi
if ! readelf -h "${THERMO_RAW_READER_DIR}/RadiantThermoRawReader" | grep -q "$elf_machine"; then
    printf 'Thermo RAW reader does not match package architecture %s.\n' "$ARCH" >&2
    exit 2
fi
mkdir -p "${radiant_bin}/thermoraw"
cp -a "${THERMO_RAW_READER_DIR}/." "${radiant_bin}/thermoraw/"

cp AlgorithmsFFLib/libAlgorithmsFFLib.so "$radiant_lib/libAlgorithmsFFLib.so"
cp ChemLib/libChemLib.so "$radiant_lib/libChemLib.so"
cp EigenLib/libEigenLib.so "$radiant_lib/libEigenLib.so"
cp FileReadersLib/libFileReadersLib.so "$radiant_lib/libFileReadersLib.so"
cp UtilsLib/libUtilsLib.so "$radiant_lib/libUtilsLib.so"
cp WorkFlowsLib/libWorkFlowsLib.so "$radiant_lib/libWorkFlowsLib.so"
cp PyTorchLib/libPyTorchLib.so "$radiant_lib/libPyTorchLib.so"

TIMSREADER_ARCH_DIR="${TIMSREADER_ARCH_DIR:-}"
if [ -z "${TIMSREADER_ARCH_DIR}" ]; then
    case "${ARCH}" in
        amd64)
            TIMSREADER_ARCH_DIR="x86_64"
            ;;
        arm64)
            TIMSREADER_ARCH_DIR="aarch64"
            ;;
    esac
fi

TIMSREADER_RUNTIME_LIBRARY="${TIMSREADER_RUNTIME_LIBRARY:-}"
if [ -z "${TIMSREADER_RUNTIME_LIBRARY}" ] && [ -f CMakeCache.txt ]; then
    TIMSREADER_RUNTIME_LIBRARY=$(sed -n 's/^TIMSREADER_LIBRARY:FILEPATH=//p' CMakeCache.txt)
fi
if [ -z "${TIMSREADER_RUNTIME_LIBRARY}" ] && [ -n "${TIMSREADER_ARCH_DIR}" ]; then
    TIMSREADER_RUNTIME_LIBRARY="/src/PythiaDIACpp/ThirdPartyLibs/timsreader/lib/${TIMSREADER_ARCH_DIR}/libtimsreader.so"
fi

if [ -n "${TIMSREADER_RUNTIME_LIBRARY}" ] && [ -e "${TIMSREADER_RUNTIME_LIBRARY}" ]; then
    cp "${TIMSREADER_RUNTIME_LIBRARY}" "${radiant_lib}/libtimsreader.so"
fi

TORCH_LIBRARY_DIR="${TORCH_LIBRARY_DIR:-/src/pytorch/build/lib}"
cp "${TORCH_LIBRARY_DIR}/libtorch.so" "$radiant_lib/"
cp "${TORCH_LIBRARY_DIR}/libtorch_cpu.so" "$radiant_lib/"
cp "${TORCH_LIBRARY_DIR}/libc10.so" "$radiant_lib/"

# Prebuilt libtorch bundles a libgomp runtime; source builds do not.
for f in "${TORCH_LIBRARY_DIR}"/libgomp*.so*; do
    [ -e "$f" ] || break
    cp "$f" "$radiant_lib/"
done

# ARM PyTorch wheels place OpenMP and Arm Compute dependencies outside torch/lib.
# Copy their runtime libraries only; Python modules and headers are not needed.
TORCH_EXTRA_LIBRARY_DIR="${TORCH_EXTRA_LIBRARY_DIR:-${TORCH_LIBRARY_DIR}/../../torch.libs}"
if [ -d "$TORCH_EXTRA_LIBRARY_DIR" ]; then
    for f in "${TORCH_EXTRA_LIBRARY_DIR}"/*.so*; do
        [ -e "$f" ] || continue
        cp -L "$f" "$radiant_lib/"
    done
fi

# Create symlinks
# Must match ${radiant_bin} above, but relative to the target!
ln -sfn radiant/RadiantDIA "${package_dir}/usr/local/bin/RadiantDIA"

mkdir -p "${package_dir}/DEBIAN/"
cp "control.${ARCH}" "${package_dir}/DEBIAN/control"
echo "Version: ${radiantdia_version}" >> "${package_dir}/DEBIAN/control"
printf 'activate-noawait ldconfig\n' > "${package_dir}/DEBIAN/triggers"

echo "Building RadiantDIA DEB…"

dpkg-deb --root-owner-group --build "$package_dir"
