# Direct Thermo RAW scan helper

`RadiantThermoRawReader` opens Thermo `.raw` files with the RawFileReader API
provided by the ThermoRawFileParser distribution. It calls the parser library's
public `ScanTrailer` and `SpectrumWriter.GetReaction` utilities and follows
`SpectrumWriter.ReadMZData` for native centroid streams, segmented centroid data,
and profile scan centroiding. Reference and exception centroid data are included,
matching the parser's default. It sends binary spectra directly to Radiant;
it never invokes the converter or creates XML/mzML. RAW loading is eager: the
helper emits the complete MS1/MS2 stream, and Radiant stores all emitted scans.
It does not provide random access or mzML-style lazy loading.

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
parser distribution for any target runtime. The selected directory is always
validated and used; automatic dependency preparation is performed only when no
parser directory is supplied. An ARM64 parser directory must contain
`radiant-thermo-dependencies.json` declaring AnyCPU assemblies. `--output`, `--dotnet`, and `--runtime` select alternate
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

When no parser directory is supplied, the preparation script downloads
compatible dependencies automatically. It downloads ThermoRawFileParser source commit
`cf548e4a3fadcc796291a01fdbd088af986a5d72`, rebuilds its assembly with
`PlatformTarget=AnyCPU`, and uses official Thermo RawFileReader/Data 8.0.42
packages from vendor commit `47eeda0d54a9ca28ace8e6548b7c01b6843207c2`.
These packages provide AnyCPU assemblies. Source and package SHA256 values
are checked before building; cached outputs live in
`build/thermoraw/dependencies/`. The published
`radiant-thermo-dependencies.json` records the dependency pins. Vendor assemblies
are used as supplied. First builds require network access, `curl`, `tar`,
`sha256sum`, and the .NET SDK.

The implementation uses ThermoRawFileParser source commit
`cf548e4a3fadcc796291a01fdbd088af986a5d72`. The dependency is the parser's
`ThermoRawFileParser.dll`, `ThermoFisher.CommonCore.Data.dll`,
`ThermoFisher.CommonCore.RawFileReader.dll`, and `log4net.dll`. Upstream source:
`https://github.com/compomics/ThermoRawFileParser`.

Automated tests use synthetic binary streams and do not require Thermo RAW data.
Real-RAW parity qualification is an external release check: it should compare
the current helper and dependency hashes with a paired mzML input and inspect
both scan-level data and downstream Radiant results. The branch does not claim
that synthetic protocol tests establish general Thermo RAW correctness.

The helper emits protocol 4. Numeric fields are little-endian:

| Part | Fields |
| --- | --- |
| Header | 8 bytes `RDRRAW01`, uint32 version `4`, uint32 count `0xffffffff` for streaming |
| Batch | uint32 payload byte count, uint32 scan count, then that many bytes of concatenated scan records |
| Scan | uint32 payload bytes `36 + 8 * points`; int32 native scan number in both scan-number fields, MS level; float32 retention time in minutes, collision energy, precursor target m/z, lower isolation offset, upper isolation offset; uint32 point count; pairs of float32 m/z and intensity |
| Footer | 8 bytes `RAWDONE1`, uint32 emitted spectrum count |

MS1 fields without precursor/activation metadata are `-1`. Only MS1 and MS2
spectra are included. Independent reader workers process disjoint native scan
ranges and publish bounded batches to one stdout serializer. Batches can arrive
out of order; Radiant sorts by native scan number and assigns the same sequential
scan numbers used by its mapped mzML reader. Older fixed-count protocol 2 and
ordered streaming protocol 3 are still accepted by Radiant. A scan-time filter
is applied after the complete RAW stream has been loaded.
Retention
times are cast to float minutes before conversion to seconds in Radiant, matching
the mapped mzML reader's operation order. Per-scan memory is bounded by a maximum
of ten million peaks. A missing footer or nonzero process exit indicates failure.
