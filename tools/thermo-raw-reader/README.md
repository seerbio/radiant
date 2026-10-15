# Direct Thermo RAW scan helper

`RadiantThermoRawReader` opens Thermo `.raw` files with the RawFileReader API
provided by the ThermoRawFileParser distribution. It calls the parser library's
public `ScanTrailer` and `SpectrumWriter.GetReaction` utilities and follows
`SpectrumWriter.ReadMZData` for native centroid streams, segmented centroid data,
and profile scan centroiding. Reference and exception centroid data are included,
matching the parser's default. It sends binary spectra directly to Radiant;
it never invokes the converter or creates XML/mzML.

Docker images and `.deb` packages include the helper in a `thermoraw/` directory
beside Radiant and use the shared .NET 8 runtime. Docker includes the runtime;
installing the `.deb` with APT automatically installs its `dotnet-runtime-8.0`
dependency. Packaged builds discover the helper automatically, with no manual
runtime installation or `--thermo-raw-reader` option.

For a source build on x86_64 or ARM64, run:

```bash
tools/thermo-raw-reader/build.sh
```

The script uses a local .NET SDK or installs SDK 8.0.426 under
`build/thermoraw/dotnet`, then publishes a framework-dependent Linux helper to
`build/bin/thermoraw/RadiantThermoRawReader`. It selects `linux-x64` on x86_64
hosts and `linux-arm64` on aarch64/arm64 hosts, including AWS EC2 Graviton. The
SDK installer uses the build host's architecture. The SDK includes the runtime
needed for source development. If the script installed the SDK in its default
private directory, expose that runtime before running Radiant:

```bash
export DOTNET_ROOT="$PWD/build/thermoraw/dotnet"
```

For a manual deployment to a machine without the SDK, install the .NET 8
runtime there; on Ubuntu, use `sudo apt-get install dotnet-runtime-8.0`.
For an isolated manual deployment without a shared runtime, optionally publish
with `--self-contained --output=build/thermoraw/standalone` and deploy that
complete directory. This bundles a private runtime; Docker and `.deb` packaging
use the default shared-runtime publish.

Pinned AnyCPU dependencies are prepared automatically on both architectures;
an unpacked parser distribution is not required. To use an existing x86_64
ThermoRawFileParser .NET 8 distribution instead:

```bash
tools/thermo-raw-reader/build.sh \
  --parser-dir=/path/to/ThermoRawFileParser
```

Set `THERMO_RAW_FILE_PARSER_DIR` instead of supplying `--parser-dir` if convenient.
ARM64 targets always use the compatible pinned AnyCPU dependencies.
`--output`, `--dotnet` and
`--runtime` allow other publish locations, SDK locations and runtime targets.
The input distribution assemblies are copied during publishing and retain
their existing licenses; binaries are build artifacts and are not committed.

To cross-publish the ARM64 reader from an x86_64 build host:

```bash
tools/thermo-raw-reader/build.sh \
  --runtime=linux-arm64 \
  --output=build/thermoraw/arm64
```

Deploy the complete output directory to `thermoraw/` beside the Graviton
instance's Radiant executable for automatic discovery. For another location,
select its executable with `RADIANT_THERMO_RAW_READER` or `--thermo-raw-reader`.
The ARM64 reader must run with an ARM64 Radiant build.

The default dependency preparation works on x86_64 and ARM64; older x64 parser
distributions cannot load on ARM64. The preparation script downloads
ThermoRawFileParser source commit
`cf548e4a3fadcc796291a01fdbd088af986a5d72`, rebuilds its assembly with
`PlatformTarget=AnyCPU`, and uses official Thermo RawFileReader/Data 8.0.42
packages from vendor commit `47eeda0d54a9ca28ace8e6548b7c01b6843207c2`.
These packages provide AnyCPU assemblies. Source and package SHA256 values
are checked before building; cached outputs live in
`build/thermoraw/arm-dependencies/`. The published
`radiant-arm64-dependencies.json` records the dependency pins. Vendor assemblies
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
