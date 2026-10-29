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

Radiant DIA is distributed as a Docker image and as Ubuntu/Debian packages from
the [Releases page](https://github.com/seerbio/radiant/releases). Docker images
are available for Intel (`x86_64`) and ARM / Apple Silicon (`aarch64`):

```shell
docker run --rm seerbio/radiant-dia:latest --help
```

The package names distinguish the architecture, for example:

```shell
apt-get install -y radiantdia_<version>+amd64.deb   # Intel
apt-get install -y radiantdia_<version>+arm64.deb   # ARM / Apple Silicon
RadiantDIA --help
```

The `RadiantDIA` executable requires a library, FASTA, configuration, and one or
more explicit data paths:

```shell
RadiantDIA [options] <fragment-library> <fasta> <radiant-config> <data> [data ...]
```

For example:

```shell
RadiantDIA library.tsv database.fasta config.radiantConfig sample.mzML \
  --output-folder radiant-results/
```

To process a batch with the same library and configuration:

```shell
RadiantDIA library.tsv database.fasta config.radiantConfig \
  sample1.mzML sample2.mzML sample3.mzML --output-folder radiant-results/
```

Data paths are processed sequentially in the order supplied, with no fixed
batch-size limit. The spectral library and library-derived candidate pairs are
loaded once for the process. Each file starts with its own reader, calibration,
mass-accuracy optimization, candidate scores, and classifier state; fitted
settings from one file are not used for the next file.

Native outputs retain the existing `<input-filename>.radiantDIA` naming and
format. Without `--output-folder`, results are written alongside each input.
With it, results are placed in the shared output folder. A batch whose inputs
would produce the same output path is rejected before any file is processed.
Use distinct input filenames or separate invocations for those inputs.

Add `--per-file-logs` to also save each input's Radiant processing diagnostics
to `<input-filename>.log` beside its native output:

```shell
RadiantDIA library.tsv database.fasta config.radiantConfig \
  sample1.mzML sample2.mzML --output-folder radiant-results/ --per-file-logs
```

This writes `sample1.mzML.log` and `sample2.mzML.log`. Each log contains its input
and shared library/configuration paths, processing diagnostics, and outcome.
Console diagnostics remain available for existing stdout/stderr redirection.
Library initialization and the batch summary stay in the combined console log;
messages written directly by native dependencies to stdout/stderr also stay
there. Log files are replaced on each invocation. Duplicate log destinations
are rejected before initialization, including those for missing or unsupported
inputs. A file whose requested log cannot be opened is skipped and reported as
failed; log write errors also cause a failed outcome. Later inputs continue.
Without this option, Radiant does not create per-file logs.

An unsupported, missing, or failed data file is reported as `FAILED`, and
processing continues with the remaining files. Every successful workflow is
reported as `COMPLETED`; the final batch summary lists completed and failed
counts. Exit status is zero only when all files succeed. Library or configuration
initialization errors stop the entire invocation. Fatal process failures, such
as termination by the operating system, cannot be recovered within the batch.
As in single-file runs, a successful workflow with no writable candidates may
produce no native result file.

The accepted input families are:

- fragment library: `.tsv`, or `.speclib`;
- protein database: `.fasta` or `.fas`;
- configuration: `.radiantConfig` (`.toml` files are also accepted);
- search data: `.mzML`, Thermo `.raw`, or a Bruker `.d` data directory.

Thermo `.raw` files are supported as search inputs.

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
