# Radiant DIA Development

## Setting up a development environment

The scripts used to buuild Radiant also support creating and working within
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

## Running in CLion IDE
1. Run the script from **Setting up a development environment**
2. Run the script from **Running a local build**
3. In CLion->File->Settings->Build,Execution,Deployment->CMake, add the following to the _Enviornment_ textbox:  **CMAKE_PREFIX_PATH=/Path/to/radiant/pytorch/build/share/cmake/Torch** 

## Running tests

Once the project has been built, run `ctest` in the build directory:

```
./test-local.sh
```


