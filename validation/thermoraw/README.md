# Direct Thermo RAW validation

Branch `thermoraw` was created from Radiant `master` at `43546479`.
Radiant reads Thermo RAW scan metadata and native centroid peaks through
ThermoRawFileParser's library utilities and Thermo RawFileReader. The reader
streams binary scan records directly to Radiant and creates no mzML or XML.

The paired input was PXD049412 HeLa 250 pg technical replicate 4:
`20231010_250pgHeLa_240k_20Th_40ms_DDM_50SPD_4`. The reference was the existing
mzML generated with native MS1/MS2 centroiding. Both searches used the full
HeLa spectral library, the same FASTA, and the saved
[configuration](config.radiantConfig), with eight threads. Input paths,
versions, settings, and hashes are recorded in [manifest.json](manifest.json)
and [native-stream-comparison.json](native-stream-comparison.json).

| Check | RAW | mzML | Comparison |
| --- | ---: | ---: | --- |
| Spectra | 30,318 | 30,318 | Exact metadata match |
| MS1 / MS2 spectra | 1,890 / 28,428 | 1,890 / 28,428 | Exact match |
| Centroid peak pairs | 47,485,476 | 47,485,476 | Exact float64 and float32 match |
| Result rows | 67,080 | 67,080 | Exact match |
| Result columns | 188 | 188 | Every column compared |
| Target precursors, precursor q-value < 0.01 | 28,131 | 28,131 | Exact match |
| Modified peptides, precursor q-value < 0.01 | 27,799 | 27,799 | Exact match |

The complete `.radiantDIA` Parquet files are **byte-for-byte identical**.
Both have SHA-256:

```text
f4cb39cbe9dddaf1fd89531a5f7fec2c1767bbb08918b69af22ec12543f36f73
```

The comparison used zero numeric tolerance and ignored no columns.
See [results-comparison.json](results-comparison.json) for every compared
field and all identification counts. The in-Radiant reader comparison also
matched every retained peak and every scan metadata field, with matching
SHA-256 digests in [spectra-comparison.json](spectra-comparison.json).

Six relevant test suites passed: command-line file validation, base reader,
mzML reader, reader dispatch, direct RAW reader, and configuration parsing.
The direct reader tests also check malformed/truncated streams, incorrect
counts, process failure, paths with spaces, uppercase `.RAW`, and failed reopen.
See [unit-tests.log](unit-tests.log).

Exact search commands are saved in [mzml-command.json](mzml-command.json) and
[raw-command.json](raw-command.json). Output files and full runtime logs remain
under `build/thermoraw-validation/`. The measured searches finished in about
5 minutes 39 seconds for mzML and 5 minutes 56 seconds for RAW; these were
overlapping runs, so those times are execution records rather than an isolated
performance benchmark.

ARM64 build instructions and dependency preparation are in
[the scan reader guide](../../tools/thermo-raw-reader/README.md). The changed
native C++ reader, reader dispatch, validation utilities, CLI parser, and
application entry point compiled successfully to AArch64 objects with actual
ARM64 Qt headers. Commands and results are saved in
[arm64-cpp/commands.json](arm64-cpp/commands.json) and
[arm64-cpp/results.json](arm64-cpp/results.json).

The ARM64 scan reader was also built from a fresh dependency cache and run
under QEMU on an existing Thermo RAW blank. All 5,845 spectra and 158,609 peak
pairs produced a binary stream byte-for-byte identical to the x86 reader;
see [arm64-stream-comparison.json](arm64-stream-comparison.json). The AnyCPU
libraries required by ARM64 also matched the full HeLa reference's 30,318
spectra and 47,485,476 peaks exactly; see
[portable-stream-comparison.json](portable-stream-comparison.json).
[arm64-dependencies.json](arm64-dependencies.json) records the verified source
and vendor package pins. A full Radiant search on a physical Graviton instance
was not performed.
