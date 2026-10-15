# Radiant DIA Development

## Setting up a development environment

The scripts used to buuild Radiant also support creating and working within
a local development environment, with support for incremental CMake/Make builds.
Currently this is only tested to work within a Ubuntu environment.

To set up your environment:

```shell
./install-build-deps-ubuntu.sh && ./get-or-build-libtorch.sh
```

On ARM64, the default PyTorch 2.0.1 setup extracts the C++ headers, CMake files,
and libraries from a pinned official ARM64 wheel, including its bundled runtime
libraries. To build that version from source instead, run
`PYTORCH_BUILD_FROM_SOURCE=1 ./get-or-build-libtorch.sh`. Other PyTorch versions
on ARM64 use the source build path. x86_64 builds use the prebuilt LibTorch
distribution.

To enable Bruker `.d` support locally, vendor `timsreader` into:

```text
ThirdPartyLibs/timsreader/include/timsreader.h
ThirdPartyLibs/timsreader/include/timsreader.hpp
ThirdPartyLibs/timsreader/lib/x86_64/libtimsreader.so
ThirdPartyLibs/timsreader/lib/aarch64/libtimsreader.so
```

The headers are shared across architectures. `build-local.sh`, Docker builds, and GitHub Actions pick the matching `lib/<arch>/libtimsreader.so` automatically.

## Running a local build

```shell
./build-local.sh
```

To include direct Thermo RAW input in a source build, build the scan reader:

```shell
tools/thermo-raw-reader/build.sh
```

On both x86_64 and ARM64, this prepares pinned AnyCPU Thermo dependencies and
uses a local .NET SDK or installs SDK 8.0.426 for the build. The reader passes
native scan data to Radiant without creating mzML. Its self-contained output,
including the .NET runtime, goes to `build/bin/thermoraw/`, which the executable
discovers automatically. Docker and `.deb` builds include this complete helper
directory; end users need no .NET installation or manual reader selection.
Alternatively, configure `-DTHERMO_RAW_READER_EXECUTABLE=/path/to/RadiantThermoRawReader`,
set `RADIANT_THERMO_RAW_READER`, or pass `--thermo-raw-reader` when running.

An existing x86_64 ThermoRawFileParser distribution can still be selected with
`--parser-dir=/path/to/ThermoRawFileParser` or `THERMO_RAW_FILE_PARSER_DIR`.
On EC2 Graviton, the script selects ARM64 automatically. It also supports
`--runtime=linux-arm64` for cross-publication from an x86_64 build host; ARM64
always uses the compatible pinned AnyCPU dependencies. See the
[reader build guide](tools/thermo-raw-reader/README.md) for details.

## Running in CLion IDE
1. Run the script from **Setting up a development environment**
2. Run the script from **Running a local build**
3. In CLion->File->Settings->Build,Execution,Deployment->CMake, add the following to the _Enviornment_ textbox:  **CMAKE_PREFIX_PATH=/Path/to/radiant/pytorch/build/share/cmake/Torch** 

## Running tests

Once the project has been built, run `ctest` in the build directory:

```
./test-local.sh
```

`MsReaderThermoRawTests` tests the native scan protocol and process failures
without requiring Thermo assemblies or a vendor RAW fixture. To additionally
compare every spectrum and peak in a real RAW/mzML pair:

```shell
RADIANT_TEST_THERMO_RAW=/path/to/sample.raw \
RADIANT_TEST_THERMO_MZML=/path/to/equivalent.mzML \
RADIANT_TEST_THERMO_REPORT=/path/to/spectra-comparison.json \
  build/bin/MsReaderThermoRawTests realRawMzMLExactEquality
```

The reference mzML must have the same native MS1/MS2 centroiding settings.
After running both inputs through Radiant with the same library, database,
configuration, and `QT_HASH_SEED=0`, compare all result columns exactly:

```shell
python3 scripts/compare_thermoraw.py mzml-results/sample.mzML.radiantDIA \
  raw-results/sample.raw.radiantDIA --output comparison.json
```

The comparison script requires Python with NumPy and PyArrow.
