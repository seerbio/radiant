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

Thermo RAW helper dependencies are built automatically by `build-local.sh` and
the Docker build. The helper build script is an implementation detail of those
builds; it does not require a separate developer or CI setup step.

## Running in CLion IDE
1. Run the script from **Setting up a development environment**
2. Run the script from **Running a local build**
3. In CLion->File->Settings->Build,Execution,Deployment->CMake, add the following to the _Enviornment_ textbox:  **CMAKE_PREFIX_PATH=/Path/to/radiant/pytorch/build/share/cmake/Torch** 

## Running tests

Once the project has been built, run `ctest` in the build directory:

```
./test-local.sh
```

To compare a shared-library batch with independent searches on real inputs,
install Python `pyarrow` and run:

```shell
python3 tools/validate-multi-file.py \
  --binary build/bin/RadiantDIA --library library.tsv --fasta database.fasta \
  --config config.radiantConfig --output validation-results \
  sample1.mzML sample2.mzML
```

Use a new output directory and inputs with distinct filenames. The helper
records commands, logs, timing, library-load counts, and full native-output
comparisons in `summary.json`. It requires identical schemas, nonempty matching
row counts, exact categorical values, and numerical agreement within absolute
and relative tolerances of `1e-6` by default. Override the tolerances with
`--absolute-tolerance` and `--relative-tolerance` when appropriate for the
validation protocol.

