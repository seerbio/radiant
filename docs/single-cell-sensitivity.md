# Single-cell sensitivity options

These options reproduce the native single-cell development policies without
environment-variable patches. They are opt-in while the full benchmark is
being completed. They do not enable match between runs.

| Section | Option | Default | Behavior |
|---|---|---:|---|
| MS1Params | alignMs1ScanTimes | false | Interpolate the MS2 reference at actual MS1 acquisition times during main non-TIMS extraction. Preserve the measured MS1 intensity trace. |
| MS2Params | mainMinSimultaneousFragments | 4 | Minimum simultaneous leading-fragment support for main non-TIMS peak extraction; accepts 3–12. |
| MS2Params | ionsSharedToReject | 4 | Minimum distinct shared supported fragments linking candidates for pre-NN competition; accepts 2–12. |
| MS2Params | postNeuralNetSharedFragments | 0 | Optional physical competition after successful non-TIMS NN scoring; zero disables it, 2–12 sets the shared-fragment threshold. |
| NeuralNetParams | neuralNetGroupPeptideFamilies | false | Keep all charges, modifications, I/L equivalents, and target/decoy members of an origin peptide in the same held-out fold. |
| NeuralNetParams | neuralNetEnsembleSize | 1 | Independently initialized networks trained inside each held-out fold. |
| NeuralNetParams | neuralNetLogIntensities | false | Apply log1p to raw intensity features before min-max scaling. Existing log-intensity features are unchanged. |
| NeuralNetParams | neuralNetShuffleEachEpoch | false | Reorder training rows and their labels together each epoch using a deterministic RNG local to each model. |

Calibration keeps the existing MS1 scoring and four-shared-fragment
competition rule. The main-extraction settings therefore do not silently
change the calibration policy. TIMS retains its existing three-fragment
seed and separate MS1 processing.

The checked-in experimental single-cell configuration uses three family
folds, four networks per fold, 24 epochs with reshuffled training rows,
a 50,000-candidate NN limit,
log-transformed intensities, corrected MS1 timing, three-fragment seeds,
three-shared-fragment pre-NN competition, and two-shared-fragment post-NN
competition. Prediction normalization is
disabled for this profile. The default configuration retains its previous
behavior apart from the empty-candidate and shorter-MS1-trace bug fixes.
MS1 calibration anchors also use the origin precursor mass for both target
and decoy candidates, matching the mass at which their MS1 traces were
extracted. This corrects the fitted anchor reference; it does not enable
the experimental application of MS1 mass corrections.

Post-NN competition uses integrated unshared fragment evidence, with the
fixed held-out NN probability breaking equal-evidence ties. It does not
change neural training, intensities, or LDA scores. After filtering, retain
one best peak for each origin/reported precursor and decoy state, group exact
probability ties, and recalculate precursor confidence from the reverse
cumulative minimum of `(decoys + 1) / targets`. Both `QValue` and
`PrecursorQValue` use this new confidence. Output uses exact NN probability
order, with decoys first for exact ties; it does not reuse the earlier
approximate NN/LDA comparator. `DecoyRatio` retains its upstream diagnostic
value and is not the final precursor q-value. Existing peptide and protein
annotations remain separate. The operation is skipped for TIMS and for
LDA fallback when NN training is unavailable.

Performance must be assessed with the same search library and no-MBR
policy. Report precursor-only and joint precursor/protein confidence
separately. Nominal 1% q-values are not an empirical guarantee of 1% false
discoveries; the development benchmark also measures entrapment error.
The DIA-NN 1.8.1 performance target has not yet been demonstrated.

Build the executable and all Radiant libraries together after changing
these settings: the parameter structure has new fields.


## Exporting candidates and scoring complementary views

`[NeuralNetParams] candidateBundleLimit = 200000` writes a new
`<raw-file-name>.radiantCandidates` directory beside the ordinary result.
The default is zero (disabled); supported nonzero limits are2..1000000.
Export supports non-TIMS data. It uses the native LDA order after pre-NN
competition and eligibility filtering, before any neural transforms. It does
not change the ordinary neural candidate limit or the ordinary report.

Each bundle contains the242 raw feature columns, origin/reported identities,
searched fragment evidence and a separate full report. Its manifest records
file hashes, the ordered feature schema, source/library/FASTA/executable
hashes and extraction settings. The reader verifies these before scoring.
Protein annotations are read only when writing the final report.

Use the `single-cell-bundle-4.radiantConfig` and
`single-cell-bundle-3.radiantConfig` profiles in separate output directories
for each raw file. Keep the four-fragment view first in the cohort's `views`
array; view order breaks equal within-view-rank ties. Then run:

```text
RadiantRescore cohort.json
```

`configs/combined-candidate-cohort.example.json` shows the input structure.
The output directory must be new. Views within one run must have the same
raw-source hash; all runs must have matching library and FASTA hashes.
The same raw source cannot be repeated under multiple run names.

For combined rescoring, add `candidateBundleOnly = true` under
`[NeuralNetParams]` to stop each search after the candidate bundle has been
written. This optional mode requires a nonzero bundle limit and does not
write an ordinary per-view `.radiantDIA` report. The exported candidate
payload and subsequent `RadiantRescore` inputs are unchanged. A bundle
without eligible candidates is reported as an error in this mode.

The cohort JSON also accepts an optional integer `threads` from 1 to 64
(default 3). It controls the number of independent neural models trained
concurrently. Seeds, training rows, family folds, model settings, and
ensemble summation order are unchanged. For twelve available CPU workers,
use `"threads": 12`; all three folds and four ensemble models can then
train concurrently.

`RadiantDIA --candidate-cohort cohort.json` can generate these bundles in one
process. See `configs/candidate-cohort.example.json`. Its base configuration
must enable `candidateBundleOnly`; each view changes only the minimum and
shared fragment counts. Paths are relative to the cohort JSON. Output
directories must be new. Parsed library rows are reused across files, and
raw loading, calibration and tolerance fitting are reused between views of
one file. Every view still performs its own full search and competition.
Each file has fresh mutable workflow state. Library preparation is part of
the command's runtime and must be included in performance comparisons.

The fixed combined policy trains three family folds with four networks each,
24 epochs, log-intensity transformation and deterministic row shuffling. It
uses100000 candidates per run for the per-run and pooled original models,
and200000 candidates per run for the expanded pooled model. On the original
100000 candidates, it first averages the logits of the original pooled and
per-run probabilities, then averages that result with the expanded pooled
logit. The other candidates retain the expanded pooled probability.
It selects one best peak per precursor and decoy state before
physical two-fragment competition and fresh exact-tie +1 confidence.
Targets use their reported modified sequence and charge. Decoy identity must
also retain the origin modified sequence: synthetic display sequences can
collide for distinct precursor hypotheses. Keeping the reported sequence as
well preserves distinct library and synthetic decoys from one origin.
Display names and all measured
evidence remain unchanged. The native identity correction has separate
regression and fixed-prediction replay validation; downstream global identity
handling and its seven-file impact must be checked before final promotion.

Reports retain original measurements and protein annotations. `QValue` and
`PrecursorQValue` contain the new per-run precursor confidence. Peptide and
protein q-values are set to1 and their best flags to0: this command does not
perform protein inference or global confidence. These require a separate
matching downstream analysis. The earlier LDA `DecoyRatio` diagnostic is reset
to its unassigned value,-1. The result manifest records the model settings,
verified input-bundle manifests/file hashes and executable/runtime provenance.
A run with no surviving precursors is recorded
with a null report path. No identities, RTs or peak intensities are transferred
between runs; the cross-run operation pools model-training rows only.

The two-cell fresh raw benchmark is complete. Both ordinary three-fragment
reports reproduce the previous V5 reports byte for byte. An independent
decoder verified all 800,000 exported candidates, and the combined reports
preserve all 194 source columns outside the ten newly assigned confidence
and classification fields. Twelve native suites/probes pass.

At nominal 1% per-run precursor confidence, the combined B18/B20 union is
9,860 human precursors and 194 paired-shuffle entrapments, compared with
9,430/178 for V5 and 10,449/195 for the matched DIA-NN 1.8.1 search. The
estimated combined FDP is 3.85%, versus 3.70% and 3.66%, respectively.
Matched downstream processing gives a confidence-only joint precursor/protein
union of 9,232/69, versus DIA-NN's 9,888/72. These are previously examined
development cells. The frozen seven-cell V8 baseline and matched downstream
analysis are complete: only B20 and C20 reach the 95% per-file ID count target.
The seven-file precursor union is 13,394 human / 450 entrapment at per-run
q <= 0.01, or 11,225 / 77 with joint precursor/protein confidence. These are
baseline results; the DIA-NN performance target remains unmet.

MS1 precursor detection is optional for final precursor eligibility. The
high MS1 correlation requirement selects MS1 calibration anchors, not final
identifications. The current learner still uses MS1 evidence as a scoring
feature; missing or weak MS1 must be evaluated separately from a hard
requirement for precursor detection. Further candidate-detector experiments
retain an independent MS2-only view alongside any optional MS1-assisted view.
