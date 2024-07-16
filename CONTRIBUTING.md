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

## Running in CLion IDE
1. Run the script from **Setting up a development environment**
2. Run the script from **Running a local build**
3. ~~In CLion->File->Settings->Build,Execution,Deployment->CMake, add the following to the _Enviornment_ textbox:  **CMAKE_PREFIX_PATH=/Path/to/PythiaDIACpp/pytorch/build/share/cmake/Torch**~~
4. In CLion->File->Settings->Build,Execution,Deployment->CMake, add the following to the _Environment_ textbox: **-G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=/snap/clion/284/bin/ninja/linux/x64/ninja -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake** Replace the path for DCMAKE_MAKE_PROGRAM with the default value shown.

## Running tests

Once the project has been built, run `ctest` in the build directory:

```
./test-local.sh
```


