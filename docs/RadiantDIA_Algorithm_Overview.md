# Radiant DIA Algorithm Overview

This document summarizes how Radiant DIA works based on the current source code in the `radiant` repository. The implementation still uses several legacy `Pythia` class and method names, but the executable exposed to users is `RadiantDIA`.

## Executive Summary

Radiant DIA is a library-centric DIA search and scoring pipeline. At a high level, it:

- reads a fragment library and builds paired target/decoy peptide candidates
- reads DIA mass-spectrometry data and organizes MS1/MS2 signal into scan-time-aware frames
- calibrates retention time and mass accuracy from an initial high-confidence subset
- optimizes extraction tolerances
- scores candidate target/decoy pairs using chromatographic, spectral, isotope, and interference features
- estimates confidence first with a linear discriminant model and then with a bagged neural-network reranker
- optionally re-annotates peptides to protein groups from FASTA
- writes a Parquet-based `.radiantDIA` result file

The main workflow entry point is `PythiaDIAFFWorkflow::processFile()` in `WorkFlowsLib/src/PythiaDIAFFWorkflow.cpp`.

## Inputs

The `RadiantDIA` CLI expects four positional inputs plus an optional output folder:

- a fragment library file (`.fragLibFF`, `.csv`, `.tsv`, or `.speclib`)
- a FASTA file
- a configuration file (`.radiantConfig`, `.pythiaConfig`, or `.toml`)
- an MS data file (`.prqFF`, `.mzML`, or Bruker `.d`)

This argument contract is implemented in `Apps/RadiantDIA/src/CommandLineParser.cpp`.

The configuration object is parsed into `PythiaParameters`, which controls charge range, peptide length range, extraction tolerances, calibration volume, feature thresholds, neural-net parameters, and reporting behavior.

## Stage 1: Library Ingestion and Candidate Construction

The initialization path begins in `PythiaDIAFFWorkflow::init()`:

1. validate the parameter object and required input files
2. read fragment-library rows with `FragLibReader::getFragLibReaderRows()`
3. initialize `TargetDecoyCandidatePairManager`
4. retrieve pointers to the resulting target/decoy candidates

`TargetDecoyCandidatePairManager` is responsible for converting library rows into searchable peptide candidates.

### Candidate generation details

For each fragment-library row, the manager:

- extracts peptide sequence with modifications and precursor charge from the library key
- filters out entries outside the configured peptide-length and charge-state ranges
- constructs a `TargetDecoyCandidatePair`
- computes a decoy mass delta from internal residue mutation logic
- chooses the decoy fragment-shift mode from configuration

After the initial load, the manager performs additional cleanup:

- remove null/empty sequences
- collapse leucine/isoleucine isomers by replacing `I` with `L` in the deduplication key
- mark cases where a generated decoy sequence would collide with an existing target sequence

In practice, this stage produces the peptide-centric search space that the downstream scoring stages will evaluate against the DIA data.

## Stage 2: MS Data Loading and Signal Organization

When `processFile()` starts, it opens the MS file through `MsReaderPointerAcc`, which abstracts over the underlying reader implementation.

From the file, Radiant collects:

- scan-number to scan-time mapping
- global scan-time bounds
- unique tandem-scan window definitions
- MS1 scan points
- DIA target frames grouped by `MzTargetKey`

`TargetDecoyCandidatePairScoretron2::init()` then builds structures that make extraction and scoring efficient:

- an `MsFrame` and `TurboXIC` for MS1
- per-window `MsFrame` objects for DIA MS2 target windows
- an averagine lookup table used for isotope-related features

Conceptually, this stage transforms raw scan data into indexable chromatographic and spectral frames that can be queried quickly during scoring.

## Stage 3: Calibration

Radiant next attempts to learn calibration from a subset of the data using `MsCalibratomaticSettertron`.

The calibration stage:

- selects up to 16 unique tandem scan windows for training
- builds temporary `TurboXIC` objects for those windows
- tranches the target/decoy candidate set for parallel processing
- scores calibration candidates using a dedicated calibration feature set from `DiscriminantScoretron::featuresCalibration()`
- uses a first-pass discriminant model and FDR counts to identify significant peptide-spectrum matches
- refines retention-time and mass-calibration statistics with `honeIRTAndMassCalibration()`

Important implementation choices:

- calibration uses `topNMS2IonsCalibration = 6`
- the minimum peak count for calibration scoring is `4.9`
- if enough confident IDs are found, the learned LDA weights are retained for later stages

The calibration phase is intended to make later extraction and scoring more accurate by tightening RT and mass expectations around the observed run.

## Stage 4: Extraction-Width Optimization

If the config does not explicitly override MS1/MS2 extraction widths, Radiant runs `OptimizeMassAccuracyPPMSettertron`.

This stage:

- evaluates a design-of-experiments style sweep of candidate MS2 ppm widths from `3` to `50`
- rescoring a subset of candidates at each ppm value
- computes a weighted FDR-based objective for each setting
- fits a second-order polynomial across the tested ppm values
- chooses the ppm value at the maximum of the fitted curve
- sets both `ms2ExtractionWidthPPM` and `ms1ExtractionWidthPPM` from that result

Optimization uses:

- `topNMS2Ions = 12`
- a minimum peak count of `3.9`
- `DiscriminantScoretron::featuresOptimization()`

This is a practical auto-tuning step: rather than using one hard-coded mass window for every run, Radiant estimates the best extraction tolerance from the current dataset.

## Stage 5: Main Candidate Scoring

The main search happens in `PythiaDIAFFWorkflow::mainAnalysis()`.

At this point the workflow has:

- a loaded candidate set
- calibrated or partially calibrated run metadata
- per-window MS frames and TurboXIC access
- an extraction-width setting

The scorer evaluates target/decoy candidates against DIA windows and produces `CandidateScores` objects. The scoring path pulls together:

- precursor and fragment chromatographic traces
- fragment co-elution and peak-shape consistency
- similarity between observed and theoretical fragment intensity patterns
- MS1 isotope support
- peak count and peak-length statistics
- interference and target-window location features

The implementation is heavily window-centric and parallelized. `TargetDecoyCandidatePairScoretron2` prepares per-window/per-target inputs and then calls candidate scoring logic in parallel.

## Stage 6: Linear Discriminant Ranking and First-Pass FDR

After raw candidate scoring, Radiant converts candidate scores into feature arrays and applies a linear discriminant analysis style ranking step through `PythiaDIAFFWorkflowSharedMethods::processBatch()`.

This step:

- builds aligned target/decoy score pairs
- selects the configured feature subset
- trains or updates a weight vector with `DiscriminantScoretron`
- applies the weights to produce discriminant scores
- sorts candidates by discriminant score
- assigns q-values with `QValueSettertron`
- summarizes counts at multiple FDR thresholds with `FDRCLassifierNeuralNet::outputFDRResults()`

Radiant uses different feature sets for different phases:

- `featuresCalibration()` for calibration
- `featuresOptimization()` for ppm optimization and the first-pass main ranking
- `featuresNeuralNetwork()` for the second-pass nonlinear reranker

The calibration and optimization feature sets focus on a smaller curated subset of chromatographic and spectral-quality features. The neural-network feature set is much larger and includes many raw and derived attributes.

## Stage 7: Neural-Network Reranking

Radiant then reranks candidates with a bagged neural-network ensemble in `applyNeuralNetClassifier()`.

The pipeline here is:

1. remove weak candidates with too few supporting MS2 fragments
2. normalize the selected neural-network features into `KarnnNNTarget` records
3. tranche the normalized candidates according to the configured bagging size
4. train an ensemble of `CandidateClassifier` models through `FDRCLassifierNeuralNet`
5. average predictions across the ensemble
6. write the resulting classifier score back onto each candidate
7. assign q-values again, this time using the neural-network score

Neural-network configuration comes directly from `PythiaParameters`:

- epochs
- bagging size
- learning rate
- nodes fraction
- focal-loss gamma
- thread count

This stage is the final reranking layer. The earlier discriminant score acts as a strong linear filter; the neural network adds a nonlinear refinement step on a broader feature set.

## Stage 8: Protein Annotation and Entrapment Reporting

If `reannotate` is enabled, Radiant reparses the FASTA and maps accepted peptide sequences back to protein groups with `SequenceSubstringSearchomatic`.

This stage:

- loads FASTA entries
- searches peptide sequences against the FASTA
- updates the `proteinGroup` field on each candidate
- applies a decoy naming transform for neural-net decoys

The code also contains logic for entrapment-style reporting when mixed proteomes are present. That reporting updates the `ResultsSummary` struct for both the discriminant and neural-network rankings.

## Stage 9: Final Filtering and Output

Before writing results, Radiant either:

- filters to the configured FDR threshold when decoys are not being reported, or
- retains a broader list with decoys when `reportDecoys` is enabled

The output writer serializes candidates to a Parquet-backed `.radiantDIA` file. The file name defaults to the input MS data path with the `.radiantDIA` suffix appended, or to the specified output folder if one was supplied.

There are two output shapes:

- full report: `CandidateScoresReaderRow`
- truncated report: `CandidateScoresReaderRowTrunc`

## What the Algorithm Is, Conceptually

Radiant is best understood as a multi-stage DIA identification and reranking system:

- it starts from a peptide/fragment library rather than a de novo spectral search
- it constructs explicit target/decoy candidate pairs for empirical error estimation
- it extracts chromatographic and spectral evidence in a window-aware manner
- it calibrates itself to the run before full scoring
- it uses a linear discriminant stage for efficient first-pass separation
- it uses a neural-network ensemble for second-pass reranking
- it emits structured Parquet results instead of a plain-text tabular report

So the core algorithm is not just “score spectra against a library.” It is:

1. candidate generation
2. chromatogram-aware signal extraction
3. calibration
4. feature engineering
5. target/decoy discriminant ranking
6. neural-network reranking
7. protein annotation and structured reporting

## Current Implementation Caveat Observed During Local Testing

During a local run with the bundled test assets, the executable started and processed input successfully but did not complete end to end.

From tracing the code and runtime behavior:

- calibration can fail with “Could not find enough significant PSMs for calibration”
- the workflow logs that warning and tries to continue
- later stages still depend on calibration-derived state or on non-empty scored candidate batches
- the sample run then aborts with “Did not run completely”

That does not change the intended algorithm above, but it means the bundled sample inputs currently behave more like internal test fixtures than a guaranteed successful end-to-end demo.

## Key Source Files

- `Apps/RadiantDIA/main.cpp`
- `Apps/RadiantDIA/src/CommandLineParser.cpp`
- `WorkFlowsLib/src/PythiaDIAFFWorkflow.cpp`
- `WorkFlowsLib/src/PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.cpp`
- `WorkFlowsLib/src/PythiaDIAFFWorkflowAlgos/OptimizeMassAccuracyPPMSettertron.cpp`
- `AlgorithmsFFLib/src/TargetDecoyCandidatePairManager.cpp`
- `AlgorithmsFFLib/src/TargetDecoyCandidatePairScoretron.cpp`
- `AlgorithmsFFLib/src/DiscriminantScoretron.cpp`
- `AlgorithmsFFLib/src/FDRCLassifierNeuralNet.cpp`

## Bottom Line

Radiant DIA is a calibrated, library-driven DIA search engine with explicit target/decoy modeling, chromatogram-aware scoring, LDA-style first-pass ranking, and neural-network reranking. The code is organized so that each stage narrows uncertainty before the next one runs: library filtering reduces the search space, calibration improves extraction accuracy, discriminant ranking finds the likely IDs, and the neural-network ensemble performs the final rerank before export.
