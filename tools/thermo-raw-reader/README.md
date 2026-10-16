# Direct Thermo RAW scan helper

`RadiantThermoRawReader` opens Thermo `.raw` files with the RawFileReader API
provided by the ThermoRawFileParser distribution. It calls the parser library's
public `ScanTrailer` and `SpectrumWriter.GetReaction` utilities and follows
`SpectrumWriter.ReadMZData` for native centroid streams, segmented centroid data,
and profile scan centroiding. Reference and exception centroid data are included,
matching the parser's default. It sends binary spectra directly to Radiant;
it never invokes the converter or creates XML/mzML.

Build the helper as part of a normal local or Docker build:

```bash
tools/thermo-raw-reader/build.sh
```

`build-local.sh` and Docker invoke the script automatically. Install the .NET 8
SDK with `install-build-deps-ubuntu.sh`; the helper then prepares pinned AnyCPU
dependencies and publishes a framework-dependent helper to
`build/bin/thermoraw/RadiantThermoRawReader`. It selects `linux-x64` on x86_64
hosts and `linux-arm64` on aarch64/arm64 hosts. The target machine needs the
shared .NET 8 runtime; `.deb` packages declare `dotnet-runtime-8.0` and Docker
installs it with the package.

Set `THERMO_RAW_FILE_PARSER_DIR` or pass `--parser-dir=DIR` to use an existing
parser distribution. `--output`, `--dotnet`, and `--runtime` select alternate
publish locations, SDKs, and target runtimes. `build/` and the helper's
`bin/obj` directories are build outputs and must not be committed.

To cross-publish the ARM64 reader from an x86_64 build host:

```bash
tools/thermo-raw-reader/build.sh \
  --runtime=linux-arm64 \
  --output=build/thermoraw/arm64
```

Deploy the complete output directory to the Graviton instance, then select
its executable with `RADIANT_THERMO_RAW_READER` or `--thermo-raw-reader`.
The ARM64 reader must run with an ARM64 Radiant build.

ARM64 builds prepare compatible dependencies automatically instead of using an
x64 distribution. The preparation script downloads ThermoRawFileParser source commit
`cf548e4a3fadcc796291a01fdbd088af986a5d72`, rebuilds its assembly with
`PlatformTarget=AnyCPU`, and uses official Thermo RawFileReader/Data 8.0.42
packages from vendor commit `47eeda0d54a9ca28ace8e6548b7c01b6843207c2`.
These packages provide AnyCPU assemblies. Source and package SHA256 values
are checked before building; cached outputs live in
`build/thermoraw/dependencies/`. The published
`radiant-thermo-dependencies.json` records the dependency pins. Vendor assemblies
are used as supplied. First builds require network access, `curl`, `tar`,
`sha256sum`, and the .NET SDK.

The ARM64 helper was run under QEMU on a real RAW file and produced the same
binary metadata and peak stream as the x86_64 helper, byte for byte. The AnyCPU
dependencies were also checked against an existing mzML containing 30,318
spectra and 47,485,476 peak pairs, with no metadata or peak differences.

The implementation was checked against ThermoRawFileParser 2.0.0.0 and upstream
source commit `cf548e4a3fadcc796291a01fdbd088af986a5d72`. The dependency is the
parser's `ThermoRawFileParser.dll`, `ThermoFisher.CommonCore.Data.dll`,
`ThermoFisher.CommonCore.RawFileReader.dll`, and `log4net.dll`. Upstream source:
`https://github.com/compomics/ThermoRawFileParser`.

Protocol 1 uses little-endian numeric fields:

| Part | Fields |
| --- | --- |
| Header | 8 bytes `RDRRAW01`, uint32 version `1`, uint32 MS1/MS2 spectrum count |
| Scan | uint32 payload bytes `36 + 16 * points`; int32 sequential scan number, native scan number, MS level; float32 retention time in minutes, collision energy, precursor target m/z, lower isolation offset, upper isolation offset; uint32 point count; pairs of float64 m/z and intensity |
| Footer | 8 bytes `RAWDONE1`, uint32 emitted spectrum count |

MS1 fields without precursor/activation metadata are `-1`. Only MS1 and MS2
spectra are included. Radiant receives the sequential scan number used by its
mapped mzML reader; the native scan number is carried for diagnostics. Retention
times are cast to float minutes before conversion to seconds in Radiant, matching
the mapped mzML reader's operation order. Per-scan memory is bounded by a maximum
of ten million peaks. A missing footer or nonzero process exit indicates failure.
