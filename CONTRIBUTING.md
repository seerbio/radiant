# Radiant DIA Development

## Setting up a development environment

The scripts used to buuild Radiant also support creating and working within
a local development environment, with support for incremental CMake/Make builds.
Currently this is only tested to work within a Ubuntu environment.

To set up your environment:

```shell
./install-build-deps-ubuntu.sh && ./get-or-build-libtorch.sh
```

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

To include direct Thermo RAW input, build the optional scan reader from a
ThermoRawFileParser .NET 8 distribution:

```shell
tools/thermo-raw-reader/build.sh --parser-dir=/path/to/ThermoRawFileParser
```

The reader uses the parser's assemblies directly and passes native scan data
to Radiant without creating mzML. Its self-contained output goes to
`build/bin/thermoraw/`, which the executable discovers automatically.
Alternatively, configure `-DTHERMO_RAW_READER_EXECUTABLE=/path/to/RadiantThermoRawReader`,
set `RADIANT_THERMO_RAW_READER`, or pass `--thermo-raw-reader` when running.

On EC2 Graviton, use `tools/thermo-raw-reader/build.sh` without `--parser-dir`.
The build selects ARM64, verifies and builds pinned AnyCPU dependencies, and
publishes a self-contained ARM64 reader. The same script supports
`--runtime=linux-arm64` for cross-publication from an x86_64 build host.

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
