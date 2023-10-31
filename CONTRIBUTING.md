# PythiaDIA Development

## Setting up a development environment

The scripts used to buuild Pythia also support creating and working within
a local development environment, with support for incremental CMake/Make builds.
Currently this is only tested to work within a Ubuntu environment.

To set up your environment:

```shell
./install-build-deps-ubuntu.sh && ./get-or-build-libtorch.sh
```

## Running a local build

```shell
./build-local.sh
```

## Running tests

Once the project has been built, run `ctest` in the build directory:

```
./test-local.sh
```
