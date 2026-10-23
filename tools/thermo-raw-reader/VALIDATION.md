# Thermo RAW qualification

The repository tests use synthetic helper streams so normal builds do not
require Thermo RAW data or proprietary dependencies. A release that enables
direct RAW input should additionally perform a real-data qualification outside
the repository's automated test suite.

Use a RAW file and a paired mzML file representing the same acquisition. Record
the input hashes, helper hash, dependency manifest, .NET runtime/SDK versions,
and target runtime. Run the same Radiant configuration against both inputs.

The qualification should compare, at minimum:

- MS1/MS2 scan counts and scan ordering;
- every scan metadata field consumed by Radiant;
- every retained m/z and intensity value, using exact comparison where the
  inputs and conversion path permit it;
- downstream result row counts, schemas, and scientific result fields.

For ARM64, repeat the helper stream comparison with the ARM64 published helper
and AnyCPU dependency set. A QEMU check can establish helper portability, but a
physical ARM64 Radiant run is the stronger release check.

The qualification record belongs with the release or validation project. Do
not commit RAW files, generated result files, figures, or compiled helper
outputs to this source tree. Synthetic protocol tests remain the appropriate
automated coverage for malformed streams, process failures, filtering, and
timeout handling.
