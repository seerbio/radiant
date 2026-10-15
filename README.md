<p align="center"><img alt="Radiant DIA logo" src="./static/radiant-logo-lockup.png" style="width: 40%;"></p>

**Radiant DIA**™ is a fast, peptide-centric search engine for bottom-up DIA-MS
data. It searches measured fragment-ion data against a spectral library and
protein database, scores candidate features, and produces native search
results for downstream analysis.

Most users should run Radiant DIA through the complete analysis pipeline. For a ready-to-use desktop
or containerized entry point, see the
[Radiant-Fulcrum GUI](https://github.com/seerbio/radiant-fulcrum-gui/) or the
[Radiant-Fulcrum container](https://github.com/seerbio/radiant-fulcrum-container/).
This repository is primarily for users integrating, developing, or diagnosing
the Radiant DIA search stage itself.

## Where Radiant DIA fits

The normal role of Radiant DIA is:

```text
combined analysis pipeline
        │
        ├── Radiant DIA search  ◄── this repository
        │
        └── downstream analysis and reporting
```

Radiant DIA produces intermediate search results rather than a complete
protein-level analysis or final downstream report. Running `RadiantDIA` directly
is useful for development, integration, and diagnosis; production analyses
will usually continue through the surrounding pipeline.

## Recommended entry points

- **Desktop users:** start with the [Radiant-Fulcrum GUI](https://github.com/seerbio/radiant-fulcrum-gui/),
  which provides a complete local workflow.
- **CLI or server users:** start with the [Radiant-Fulcrum container](https://github.com/seerbio/radiant-fulcrum-container/),
  which provides a reproducible containerized workflow.
- **Radiant developers and integrators:** continue below for the native
  executable, supported inputs, and configuration basics.

## Native executable

Radiant DIA is distributed as a Docker image and as Ubuntu `.deb` packages from
the [Releases page](https://github.com/seerbio/radiant/releases). Docker images
are available for Intel (`x86_64`) and ARM / Apple Silicon (`aarch64`):

```shell
docker run --rm seerbio/radiant-dia:latest --help
```

The package names distinguish the architecture, for example:

```shell
apt-get install -y ./radiantdia_<version>+amd64.deb   # Intel
apt-get install -y ./radiantdia_<version>+arm64.deb   # ARM64 Linux
RadiantDIA --help
```

The package build in this repository targets Ubuntu 22.04. APT resolves its
dependencies, including the shared .NET 8 runtime, from Ubuntu's repositories.

The `RadiantDIA` executable requires four positional arguments:

```shell
RadiantDIA [options] <fragment-library> <fasta> <radiant-config> <data>
```

For example:

```shell
RadiantDIA library.tsv database.fasta config.radiantConfig sample.mzML \
  --output-folder radiant-results/
```

The accepted input families are:

- fragment library: `.tsv`, or `.speclib`;
- protein database: `.fasta` or `.fas`;
- configuration: `.radiantConfig` (`.toml` files are also accepted);
- search data: `.mzML`, Thermo `.raw`, or a Bruker `.d` data directory.

Thermo RAW input uses the direct scan reader built from the
ThermoRawFileParser library and its Thermo RawFileReader dependencies. It reads
native centroid spectra and scan metadata into Radiant without converting to
mzML. Docker images and `.deb` packages include the reader beside Radiant and
use the shared .NET 8 runtime. Installing the `.deb` with APT automatically
installs its `dotnet-runtime-8.0` dependency; Docker includes the runtime.
Packaged builds discover the reader automatically, with no manual runtime
installation or `--thermo-raw-reader` option:

```shell
RadiantDIA library.tsv database.fasta config.radiantConfig sample.raw \
  --output-folder radiant-results/
```

Source builders run `tools/thermo-raw-reader/build.sh` once. On both x86_64 and
ARM64, it prepares pinned AnyCPU Thermo dependencies and publishes the reader
to `build/bin/thermoraw/`, where Radiant finds it automatically. Source builds
use the runtime provided by the .NET SDK; manual deployments need a .NET 8
runtime on the target host. For another reader location, pass
`--thermo-raw-reader /path/to/RadiantThermoRawReader` or set
`RADIANT_THERMO_RAW_READER`.
RAW input reads MS1/MS2 scans with native centroiding and uses eager loading;
`useLazyLoading` applies to mzML input. See the
[reader build instructions](tools/thermo-raw-reader/README.md) for dependencies
and the binary scan interface.

On EC2 Graviton, `tools/thermo-raw-reader/build.sh` selects Linux ARM64
automatically. Build Radiant natively on the ARM64 host using the normal
development instructions. To cross-publish just the scan reader from x86_64,
use `--runtime=linux-arm64`.

Run `RadiantDIA --help` for the options available in the installed build. When
running in a container, paths must be valid in the execution environment—for
example, they must refer to mounted container paths rather than host path
names.

## Configuration

Radiant configuration files are TOML-style files. The repository includes an
[example configuration](example_params.radiantConfig). Use it as a starting
point. It's recommended to keep custom configuration files with their runs
so that searches can be reproduced.

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for build and contribution instructions.
The build also vendors architecture-specific `timsreader` libraries:

```text
ThirdPartyLibs/timsreader/include/timsreader.h
ThirdPartyLibs/timsreader/include/timsreader.hpp
ThirdPartyLibs/timsreader/lib/x86_64/libtimsreader.so
ThirdPartyLibs/timsreader/lib/aarch64/libtimsreader.so
```
