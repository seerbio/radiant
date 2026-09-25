# Bruker Branch Code Changes

This document summarizes the code changes currently present on the `bruker`
branch. The central theme is adding Bruker TIMS ion mobility support to
Radiant while preserving the existing non-ion-mobility DIA behavior.

## High-Level Summary

- Added Bruker TIMS data support through the reader stack, parquet conversion,
  and downstream scoring path.
- Added AlphaDIA-inspired ion mobility processing for TIMS DIA data:
  - library-centered broad mobility extraction,
  - local RT x mobility evidence selection,
  - tighter target mobility extraction,
  - MS2 mobility profile features,
  - classifier feature integration when ion mobility is available.
- Added candidate-level precursor, peptide, and protein q-value annotation.
- Expanded spectral-library TSV parsing for ion mobility, RT, and protein group
  aliases.
- Updated the parquet schema and conversion workflow so Bruker `.d` data can be
  converted and re-read with TIMS metadata preserved.
- Added tests for TIMS ion mobility indexing and id-level q-value annotation.
- Kept IMS processing guarded so files without ion mobility do not receive IMS
  scoring changes.

## Build System Changes

### `AlgorithmsFFLib/CMakeLists.txt`

- Added the new `IdLevelQValueAnnotator` source and header.
- Added the new `TimsMs2IonMobilityIndex` source and header.

### `AlgorithmsFFLib/tests/CMakeLists.txt`

- Added `IdLevelQValueAnnotatorTests`.
- Added `TimsMs2IonMobilityIndexTests`.

### `FileReadersLib/CMakeLists.txt`

- Enabled Bruker TIMS reader compilation.
- Added TIMS SDK include paths and library linkage.
- Added SQLite include paths needed by the TIMS reader.
- Added RPATH settings so the TIMS runtime library can be found from built
  binaries.

## Candidate Scoring Changes

### `AlgorithmsFFLib/src/CandidateScorertron.h`

- Extended `init(...)` to accept:
  - `MsReaderPointerAcc *`
  - `TimsMs2IonMobilityIndex *`
- Stored those pointers on the scorer so scoring can access frame RT and TIMS
  MS2 mobility evidence.
- Passed the target candidate into matrix/XIC initialization so library ion
  mobility can influence extraction only when available.
- Added IMS-specific scoring helpers:
  - `setLibraryIonMobilityRelatedScores(...)`
  - `setMs2IonMobilityRelatedScores(...)`

### `AlgorithmsFFLib/src/CandidateScorertron.cpp`

- Added guarded TIMS ion mobility extraction. The scorer only enters the IMS
  path when:
  - a target candidate is present,
  - the target library entry has positive ion mobility,
  - a `TimsMs2IonMobilityIndex` is initialized,
  - TIMS MS2 evidence is available.
- Preserved the existing `getXICs(...)` behavior for non-IMS data.
- Added AlphaDIA-inspired constants:
  - broad mobility tolerance: `0.1` 1/K0,
  - target mobility tolerance: `0.04` 1/K0,
  - mobility FWHM smoothing scale: `0.01` 1/K0,
  - RT FWHM smoothing scale: `5.0` seconds.
- Added symmetric-profile boundary logic similar to AlphaDIA's profile-limit
  selection.
- Added local RT x mobility evidence selection for TIMS MS2:
  - extracts fragment evidence in a broad library-centered mobility window,
  - aggregates evidence by frame and ion mobility index,
  - converts frame index to RT and mobility index to drift time,
  - applies Gaussian-style RT and mobility weighting,
  - selects the strongest local mobility center,
  - extracts final MS2 XICs from a tighter mobility window around that local
    center.
- Added MS1 library mobility scoring:
  - searches TIMS MS1 frame data near precursor m/z and library mobility,
  - records the observed mobility index and drift time,
  - records library/found mobility deltas,
  - records `Ms1IntensityFoundApex100IM`.
- Added MS2 mobility profile scoring:
  - builds fragment-level mobility profiles across the candidate RT window,
  - computes weighted observed mobility and delta features,
  - computes fragment apex mobility delta statistics,
  - computes matched ion fraction,
  - computes FWHM-like mobility profile statistics over the observed local
    profile window, including empty bins where possible.
- Allows MS2 mobility evidence to replace the MS1-derived observed mobility
  when stronger MS2 evidence is available.
- Populates `proteinGroup` on candidate scores from the target/decoy pair.
- Added a guard for empty candidate score-pair containers before downstream
  scoring.

## Candidate Score Schema Changes

### `AlgorithmsFFLib/src/CandidateScores.h`

- Added IMS feature enum entries used by calibration, neural-network scoring,
  and output serialization:
  - `IonMobilityDelta`
  - `IonMobilityDeltaAbs`
  - `IonMobilityPdAbs`
  - `Ms2IonMobilityWeightedDelta`
  - `Ms2IonMobilityWeightedDeltaAbs`
  - `Ms2IonMobilityApexDeltaAbsMean`
  - `Ms2IonMobilityApexDeltaAbsStDev`
  - `Ms2IonMobilityMatchedIonFraction`
  - `Ms2IonMobilityFwhmMean`
  - `Ms2IonMobilityFwhmStDev`
- Added candidate-level q-value fields:
  - `precursorQValue`
  - `peptideQValue`
  - `proteinQValue`
- Added candidate-level best-row flags:
  - `isBestPrecursorCandidate`
  - `isBestPeptideCandidate`
  - `isBestProteinCandidate`
- Added output columns for:
  - library mobility,
  - found mobility,
  - mobility index and extraction bounds,
  - IMS delta/profile features,
  - q-values at precursor, peptide, and protein level,
  - best-candidate flags.
- Updated serialization and deserialization so the new fields round-trip through
  the candidate score table.

## Classifier Changes

### `AlgorithmsFFLib/src/DiscriminantScoretron.cpp`

- Added IMS features to the score calibration and optimization feature sets.
- Added IMS features to the neural-network input feature list.
- Added default feature weights so IMS evidence can influence classifier score
  when available:
  - rewards stronger MS2 mobility matched-ion fraction,
  - penalizes larger mobility deltas,
  - penalizes broader or more variable mobility profiles.

IMS features remain inactive for non-IMS data because their values stay at the
existing default/inactive values when no library mobility or TIMS index is
available.

## Target/Decoy Scoring Changes

### `AlgorithmsFFLib/src/TargetDecoyCandidatePairScoretron.cpp`

- Detects whether the scoring context has usable library ion mobility and TIMS
  MS2 frame data.
- Builds a `TimsMs2IonMobilityIndex` for each usable target key when TIMS data
  is available.
- Passes the reader pointer and TIMS index into `CandidateScorertron`.
- Leaves non-IMS candidate scoring on the original path.

### `AlgorithmsFFLib/tests/TargetDecoyCandidatePairTests.cpp`

- Updated expectations to account for protein group propagation on candidate
  scores.

## New TIMS MS2 Ion Mobility Index

### `AlgorithmsFFLib/src/TimsMs2IonMobilityIndex.h`

### `AlgorithmsFFLib/src/TimsMs2IonMobilityIndex.cpp`

- Added a reusable index over TIMS MS2 points.
- Initializes from:
  - target-keyed TIMS MS2 frames,
  - MS frame metadata,
  - the full ion mobility index to drift-time map.
- Stores sorted indexed points containing:
  - m/z,
  - intensity,
  - frame index,
  - ion mobility index,
  - drift time.
- Keeps the full scan-index to drift-time map, including bins with no observed
  MS2 points. This is needed for AlphaDIA-like profile windows and FWHM-style
  summaries.
- Exposes:
  - initialization status,
  - point count,
  - drift-time lookup from ion mobility index,
  - bounded extraction by m/z, frame, and mobility index.

### `AlgorithmsFFLib/tests/TimsMs2IonMobilityIndexTests.cpp`

- Verifies initialization and point counting.
- Verifies extraction by m/z, frame, and ion mobility bounds.
- Verifies drift-time lookup for an unobserved mobility bin retained from the
  full scan-index map.

## Id-Level Q-Value Annotation

### `AlgorithmsFFLib/src/IdLevelQValueAnnotator.h`

### `AlgorithmsFFLib/src/IdLevelQValueAnnotator.cpp`

- Added post-classifier q-value annotation at three identification levels:
  - precursor,
  - peptide,
  - protein.
- Collapses target and decoy candidates by level-specific keys:
  - precursor key: modified peptide plus charge,
  - peptide key: modified peptide,
  - protein key: cleaned protein group.
- Uses classifier score when available, falling back to discriminant score when
  classifier score is unavailable.
- Annotates only target candidates with reported q-values.
- Adds best-candidate flags for target and decoy candidates within each
  collapsed identification level.

### `AlgorithmsFFLib/tests/IdLevelQValueAnnotatorTests.cpp`

- Verifies precursor-, peptide-, and protein-level q-value annotation.
- Verifies best-candidate flag assignment.

### `WorkFlowsLib/src/PythiaDIAFFWorkflow.cpp`

- Runs `IdLevelQValueAnnotator::annotate(...)` before decoy filtering and output
  writing so result files include identification-level q-values and best-row
  flags.

## Bruker TIMS Reader Changes

### `FileReadersLib/src/MsReaderBase.h`

### `FileReadersLib/src/MsReaderBase.cpp`

- Added TIMS MS2 frame data structures.
- Added native Bruker frame and scan number fields to scan metadata.
- Added accessors for:
  - target-keyed TIMS MS2 frame maps,
  - ion mobility index to drift-time maps.
- Ensured TIMS-specific data structures are cleared on reader reset/close.

### `FileReadersLib/src/MsReaderBrukerTims.cpp`

- Added support for reading Bruker DIA frame/window metadata.
- Reads frame-to-window-group mappings from Bruker metadata and falls back to
  cyclic DIA window grouping when needed.
- Builds TIMS MS1 and MS2 maps with frame, scan, ion mobility index, drift time,
  m/z, and intensity information.
- Preserves native Bruker frame and scan identifiers in scan metadata.
- Uses the full `scanNumToOneOverK0` map for mobility lookup.
- Trims edge scans from DIA isolation windows when possible to reduce boundary
  artifacts.
- Builds per-target-key TIMS MS2 frame maps that downstream IMS scoring can use.

## Parquet Reader and Converter Changes

### `FileReadersLib/src/MsReaderParquet.h`

### `FileReadersLib/src/MsReaderParquet.cpp`

- Extended parquet read/write support for TIMS data.
- Preserves and reconstructs:
  - ion mobility drift time,
  - ion mobility index,
  - native frame number,
  - native scan number,
  - TIMS MS1 frame maps,
  - TIMS MS2 target-keyed frame maps.
- Writes per-mobility-bin rows for TIMS data so `.prqFF` files retain enough
  information for downstream IMS scoring.

### `FileReadersLib/src/ParquetReader.cpp`

- Updated parquet row handling for the expanded schemas used by TIMS-aware
  reader output and candidate score output.

### `FileReadersLib/tests/MsReaderParquetTests.cpp`

- Added coverage for ion mobility and native frame/scan metadata round-tripping.

### `WorkFlowsLib/src/ConvertMzMLToParquetWorkFlow.h`

### `WorkFlowsLib/src/ConvertMzMLToParquetWorkFlow.cpp`

- Renamed the workflow input concept from mzML-only to generic MS data input.
- Switched conversion to use `MsReaderPointerAcc` rather than directly creating
  an mzML reader.
- Allows Bruker `.d` directories to be converted to `.prqFF`.
- Produces output as `<absolute input path>.prqFF`.

### `Apps/LycaonMzMLTransmorgitron/main.cpp`

### `Apps/LycaonMzMLTransmorgitron/src/CommandLineParser.h`

### `Apps/LycaonMzMLTransmorgitron/src/CommandLineParser.cpp`

- Renamed the positional argument from `mzML-path` to `ms-data-path`.
- Accepts either mzML files or Bruker `.d` directories.
- Passes the generic MS data path through to the conversion workflow.

## Spectral Library Reader Changes

### `FileReadersLib/src/FragLibReaderRow.h`

- Added `proteinGroups` to spectral library row data.

### `FileReadersLib/src/FragLibTsvReader.h`

### `FileReadersLib/src/FragLibTsvReader.cpp`

- Added retention-time aliases, including:
  - `AverageExperimentalRetentionTime`
- Added ion mobility aliases, including:
  - `PrecursorIonMobility`
  - `IMS time`
  - `IMSTime`
  - `IM`
  - `Mobility`
- Added protein group aliases, including:
  - `ProteinGroup`
  - `ProteinId`
  - `ProteinName`
  - `UniprotID`
  - `Genes`
  - `GeneName`
- Propagates parsed ion mobility and protein group values into grouped library
  rows used for candidate generation and reporting.

### `FileReadersLib/tests/FragLibTsvReaderTests.cpp`

- Added tests for ion mobility aliases.
- Added tests for the `AverageExperimentalRetentionTime` RT alias.
- Updated assertions for grouped library rows carrying mobility information.

## Reader Pointer Accumulation

### `FileReadersLib/src/MsReaderPointerAcc.cpp`

- Updated reader selection/initialization so the common reader path can expose
  TIMS-aware data to workflows and scoring.
- Supports the Bruker `.d` and parquet TIMS data paths needed by the new IMS
  scoring code.

## Validation Performed

### Build and Unit Tests

```bash
cmake --build /home/ubuntu/git/radiant/build --target TimsMs2IonMobilityIndexTests RadiantDIA -j 16
/home/ubuntu/git/radiant/build/bin/TimsMs2IonMobilityIndexTests
```

Result:

- `TimsMs2IonMobilityIndexTests`: 3 passed, 0 failed.
- `RadiantDIA` built successfully.

### IMS Small-Window Run

Input data:

- Library:
  `/home/ubuntu/se_bakeoff/libraries/uniprot_homosapiens_isoforms_con_20221209.predicted.speclib`
- FASTA:
  `/home/ubuntu/se_bakeoff/fasta/uniprot_homosapiens_isoforms_con_20221209.fasta`
- Bruker data:
  `/home/ubuntu/pegasus/alphadia_synchropasef_test/20231218_TIMS03_PaSk_SA_K562_syPASEF_200ng_var_IM0713_S1-E7_1_41539.d`

Output:

- `/home/ubuntu/pegasus/alphadia_synchropasef_test/radiant_ims_alphadia_profile_scoring_v4/20231218_TIMS03_PaSk_SA_K562_syPASEF_200ng_var_IM0713_S1-E7_1_41539.rt400_700.prqFF.radiantDIA`

Runtime:

- `5:07.78`
- Max RSS: `28405752 kB`

Q < 0.01 counts for the Radiant IMS run:

| Metric | Value |
| --- | ---: |
| Rows | 56626 |
| Target rows | 41832 |
| PSM rows | 21311 |
| Precursors | 1285 |
| Modified peptides | 1247 |
| Peptide sequences | 1202 |
| Proteins | 888 |
| Best-peptide mean absolute IM delta | 0.015513 |
| Best-peptide mean IM FWHM | 0.000671 |
| Best-peptide median IM FWHM | 0.000671 |
| Best-peptide nonzero IM FWHM fraction | 1.0 |

AlphaDIA comparison on the same small-window data:

| Metric | Value |
| --- | ---: |
| Rows | 1338 |
| Target rows | 1338 |
| PSM rows | 2445 |
| Precursors | 1318 |
| Modified peptides | 1246 |
| Peptide sequences | 1202 |
| Proteins | 804 |
| Best-peptide mean absolute IM delta | 0.009902 |
| Best-peptide mean IM FWHM | 0.005047 |
| Best-peptide median IM FWHM | 0.004775 |
| Best-peptide nonzero IM FWHM fraction | 1.0 |

Peptide sequence overlap at Q < 0.01:

| Set | Count |
| --- | ---: |
| Radiant peptide sequences | 1202 |
| AlphaDIA peptide sequences | 1202 |
| Shared peptide sequences | 689 |
| Radiant-only peptide sequences | 513 |
| AlphaDIA-only peptide sequences | 513 |

Generated comparison artifacts:

- `/home/ubuntu/pegasus/alphadia_synchropasef_test/radiant_alphadia_profile_scoring_v4_counts_q001.csv`
- `/home/ubuntu/pegasus/alphadia_synchropasef_test/radiant_alphadia_profile_scoring_v4_vs_alphadia_peptide_sequence_overlap_q001.csv`
- `/home/ubuntu/pegasus/alphadia_synchropasef_test/plots/radiant_alphadia_profile_scoring_v4_counts_q001.png`
- `/home/ubuntu/pegasus/alphadia_synchropasef_test/plots/radiant_alphadia_profile_scoring_v4_vs_alphadia_peptide_sequence_venn_q001.png`

### Non-IMS Guard Run

Output:

- `/home/ubuntu/pegasus/alphadia_synchropasef_test/radiant_no_ims_guard_after_profile_scoring`

Runtime:

- `3:01.76`
- Max RSS: `28250020 kB`

Q < 0.01 counts:

| Metric | Value |
| --- | ---: |
| Rows | 58578 |
| Target rows | 42944 |
| PSM rows | 19827 |
| Precursors | 1148 |
| Peptides | 1114 |
| Proteins | 752 |

IMS guard checks:

| Field | Observed behavior |
| --- | --- |
| `IonMobilityLibrary` | min/median/max `0.0`; nonzero fraction `0.0` |
| `IonMobilityFound` | min/median/max `-1.0`; nonzero fraction `0.0` |
| `IonMobilityIndex` | min/median/max `-1`; nonzero fraction `0.0` |
| `IonMobilityIndexStart` | min/median/max `-1`; nonzero fraction `0.0` |
| `IonMobilityIndexEnd` | min/median/max `-1`; nonzero fraction `0.0` |
| `IonMobilityDeltaAbs` | all `0` |
| `Ms2IonMobilityFwhmMean` | all `0` |
| `Ms2IonMobilityMatchedIonFraction` | all `0` |

This confirms the IMS code path does not activate for files/libraries without
ion mobility.

## Current Relationship to AlphaDIA

The IMS processing is now intentionally AlphaDIA-inspired, but Radiant is not a
drop-in reimplementation of AlphaDIA scoring.

Borrowed/approximated behavior:

- broad library-centered mobility search,
- local RT x mobility evidence selection,
- tighter extraction around the selected local mobility center,
- mobility profile statistics,
- mobility delta features,
- classifier use of IMS features when available.

Still Radiant-native:

- candidate generation,
- XIC matrix construction,
- non-IMS feature set,
- classifier architecture and calibration,
- q-value annotation and output schema,
- target/decoy handling,
- final reporting semantics.

Known remaining differences:

- PSM row counts are not directly comparable to AlphaDIA because the tools emit
  different row semantics and collapse/report at different stages.
- The small-window run matched AlphaDIA exactly on the number of peptide
  sequences at Q < 0.01, but the overlap was partial.
- The Radiant mobility FWHM summary remains on a different scale than AlphaDIA
  in the current implementation, even though the feature is nonzero and derived
  from a local mobility profile.

## Public Plasma TIMS IMS Validation

Input files:

- Raw: `tests/ms_data/P4777_101_15_3_R1_1_1_10579.d`
- Library: `tests/libraries/hppp-lib.predicted.speclib`
- FASTA: `tests/fasta/HPPP_Fasta.fasta`

Two IMS issues were found while comparing against DIA-NN:

1. `MsCalibratomatic` built an IM mapper, but `CandidateScorertron` was still
   using raw library iIM as the TIMS extraction and IM residual center. The
   scorer now uses calibrated IM when the mapper is fitted.
2. Bruker's C++ `TimsData` wrapper defaults to analysis-global pressure
   compensation, while AlphaRaw/DIA-NN report the no-pressure-compensation 1/K0
   axis. Radiant now opens Bruker `.d` files with `NoPressureCompensation`.

The IM residual features are kept neutral before IM calibration is fitted while
still preserving `imDriftTime` for fitting the IM mapper. Once the mapper is
available, extraction and residual scoring use the calibrated mobility center.

Q < 0.01 results on the public plasma TIMS run:

| Run | Output folder | PSMs | Precursors | Modified peptides | Proteins | IM axis check |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| DIA-NN | `tests/search_results/public_plasma_compare_speclib_rerun/diann_tims` | 2469 | 2469 | 2220 | 380 protein groups | Reference |
| Radiant 10 ppm, pressure-compensated, raw IM center | `tests/search_results/public_plasma_compare_speclib_ppm10_fullnn/radiant_tims` | 1717 | 1676 | 1509 | 293 | `IonMobilityFound - AlphaRawAxis[index]` median `+0.0216` |
| Radiant 10 ppm, pressure-compensated, calibrated IM center | `tests/search_results/public_plasma_compare_speclib_imcal_ppm10_fullnn/radiant_tims` | 1763 | 1718 | 1555 | 280 | `IonMobilityFound - AlphaRawAxis[index]` median `+0.0216` |
| Radiant 10 ppm, no-pressure, calibrated IM center | `tests/search_results/public_plasma_compare_speclib_imcal_nopressure_ppm10_fullnn/radiant_tims` | 1522 | 1486 | 1343 | 267 | `IonMobilityFound - AlphaRawAxis[index]` median `0.0000` |
| Radiant 10 ppm, no-pressure, calibrated IM center, neutral pre-cal IM residuals | `tests/search_results/public_plasma_compare_speclib_imcal_nopressure_neutralpreim_ppm10_fullnn/radiant_tims` | 1658 | 1616 | 1450 | 277 | `IonMobilityFound - AlphaRawAxis[index]` median `0.0000` |

The pressure-compensated runs produce more IDs, but their reported IM coordinate
is shifted by about 18 TIMS scans relative to AlphaRaw/DIA-NN. The no-pressure
runs are the physically aligned IMS path: on shared DIA-NN/Radiant precursors,
`Radiant IonMobilityFound - DIA-NN IM` is centered at approximately zero.

Current interpretation: IMS is now being used and reported on the AlphaRaw/DIA-NN
1/K0 axis, with calibrated-IM extraction/scoring when available. The remaining
count gap is no longer explained by an IM coordinate mismatch; it is in Radiant's
candidate ranking/classifier behavior and broader scoring differences.

### Public Plasma Drop-out Diagnostics

To verify whether the DIA-NN gap was still IMS-related, a default-off diagnostic
flag, `writeFullCandidateDebug`, was added under `[General]`. When enabled,
Radiant writes compact target-candidate TSVs before final report filtering:

- `*.radiantDIA.scored_candidates.tsv`: all scored target candidates before NN
  filtering.
- `*.radiantDIA.nn_candidates.tsv`: target candidates retained for NN scoring,
  before the q <= 0.5 report truncation.

The diagnostic run was:

- `tests/search_results/public_plasma_compare_speclib_imcal_nopressure_neutralpreim_ppm10_fullnn_debug/radiant_tims`

DIA-NN had 2469 unique precursor keys at Q < 0.01. Their farthest Radiant stage
was:

| Farthest Radiant stage | DIA-NN precursor keys |
| --- | ---: |
| Radiant final Q < 0.01 | 1480 |
| Final report, but not Q < 0.01 | 476 |
| NN-scored, but not final report | 182 |
| Scored before NN, but not retained for NN | 110 |
| Not scored by Radiant | 221 |

The DIA-NN-only keys that Radiant scored generally had lower DIA-NN evidence and
lower Radiant MS2 evidence than the shared Q < 0.01 set, but their IM residuals
were not the dominant outlier:

| Category | DIA-NN median CScore | DIA-NN median evidence | Radiant median `IonMobilityDeltaAbs` | Radiant median `Ms2IonMobilityMatchedIonFraction` |
| --- | ---: | ---: | ---: | ---: |
| Radiant final Q < 0.01 | 0.999846 | 4.258693 | 0.010604 | 1.000000 |
| Final report, not Q < 0.01 | 0.988681 | 3.192837 | 0.010907 | 0.833333 |
| NN-scored, not final | 0.985327 | 3.167624 | 0.014779 | 0.833333 |
| Scored, not NN-retained | 0.963230 | 2.997108 | 0.010989 | 0.583333 |
| Not scored | 0.927645 | 2.521482 | n/a | n/a |

At Radiant target rank 2469, there were already 338 decoys ahead of or tied into
the classifier-sorted list, so the DIA-NN count cannot be recovered by changing
only the report threshold. The remaining issue is classifier/target-decoy
separation for lower-evidence candidates.

### Public Plasma Classifier Variants

Two default-preserving neural-net config controls were added:

- `[NeuralNetParams] neuralNetCandidateLimit`, default `50000`, replacing the
  previous hard-coded NN candidate cap.
- `[NeuralNetParams] normalizeNeuralNetPredictions`, default `true`, allowing
  raw sigmoid decoy probabilities to be used directly for ranking.

Q < 0.01 variant results:

| Variant | Output folder | PSMs | Precursors | Peptides | Proteins | DIA-NN precursor overlap |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Baseline debug, normalized NN, 50k | `public_plasma_compare_speclib_imcal_nopressure_neutralpreim_ppm10_fullnn_debug/radiant_tims` | 1633 | 1591 | 1453 | 281 | 1480 |
| Normalized NN, 100k | `public_plasma_compare_speclib_variant_nn100k/radiant_tims` | 1732 | 1688 | 1528 | 291 | 1562 |
| Normalized NN, 200k | `public_plasma_compare_speclib_variant_nn200k/radiant_tims` | 1626 | 1587 | 1449 | 274 | 1473 |
| Raw NN, 50k | `public_plasma_compare_speclib_variant_rawnn50k/radiant_tims` | 1832 | 1784 | 1610 | 297 | 1641 |
| Raw NN, 100k | `public_plasma_compare_speclib_variant_rawnn100k/radiant_tims` | 1757 | 1711 | 1548 | 295 | 1591 |
| Raw NN, 50k, focal gamma 2 | `public_plasma_compare_speclib_variant_rawnn50k_focal2/radiant_tims` | 1777 | 1733 | 1562 | 292 | 1602 |

Disabling fold-wise NN prediction normalization was the best tested scoring
change on this run. Increasing the candidate cap to 100k helped under normalized
scoring, but did not combine well with raw NN scores. A 200k cap degraded the
classifier, likely by adding too many lower-quality candidates/decoys.

Current blocker: even the best tested Radiant variant remains below DIA-NN
(`1832` vs `2469` PSMs; `1784` vs `2469` precursors). Since IMS alignment and
calibrated IM features are behaving consistently, further convergence likely
requires deeper classifier/FDR work rather than additional IM-coordinate fixes.

## Review Checklist

- Confirm TIMS SDK paths are portable enough for the intended deployment
  environment.
- Review the IMS feature weights in `DiscriminantScoretron.cpp` against larger
  validation data.
- Decide whether raw NN prediction scores should become the default, or remain a
  config-controlled option pending multi-file validation.
- Decide whether `neuralNetCandidateLimit` should be surfaced in production
  configs and what dataset-scale default should be used.
- Decide whether output row semantics should be further harmonized with
  AlphaDIA for easier direct comparison.
- Add broader integration coverage for `.d -> .prqFF -> RadiantDIA` workflows
  if Bruker TIMS support is promoted beyond this branch.
