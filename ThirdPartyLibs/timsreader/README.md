# Vendored `timsreader`

Place the shared headers once and vendor one Linux shared library per CPU architecture:

```text
ThirdPartyLibs/timsreader/include/timsreader.h
ThirdPartyLibs/timsreader/include/timsreader.hpp
ThirdPartyLibs/timsreader/lib/x86_64/libtimsreader.so
ThirdPartyLibs/timsreader/lib/aarch64/libtimsreader.so
```

Notes:

- `timsreader.hpp` is required because Radiant includes the C++ wrapper API.
- `build-local.sh`, `Dockerfile`, and the GitHub Actions container jobs all select the library by the build machine architecture automatically.
- If you need to test a non-vendored checkout locally, override `TIMSREADER_ROOT`, `TIMSREADER_INCLUDE_DIR`, or `TIMSREADER_LIBRARY` in CMake.
