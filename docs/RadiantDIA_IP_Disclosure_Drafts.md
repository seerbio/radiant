# RadiantDIA IP Disclosure Drafts

This document is a drafting aid for counsel and internal review. It is not a legal opinion on patentability, freedom to operate, or validity.

## Prior-Art Framing

The current prior-art boundary for this space is crowded around:

- generic library-based DIA search
- chromatogram-library workflows
- generic target-decoy FDR
- generic neural-network reranking
- generic compute-oriented mass-spectrometry file formats
- generic audit logging and provenance capture

For that reason, the two most plausible invention surfaces are narrower:

1. An adaptive calibration-rescue controller that branches the workflow based on calibration sufficiency and emits multi-axis confidence outputs rather than treating the run as only success or failure.
2. An interference-conditioned local null model for DIA FDR that matches decoys or null cohorts to the local interference regime instead of relying only on a global null.

The point of these drafts is to claim orchestration, confidence propagation, and local error modeling rather than re-claiming broad DIA search concepts that are likely already occupied.

## Draft 1: Adaptive Calibration-Rescue Controller for DIA Analysis

### Working Title

Adaptive calibration-rescue controller for data-independent acquisition proteomics workflows.

### Problem Statement

Current DIA pipelines generally assume that enough high-confidence evidence exists to complete calibration, mass-accuracy tuning, and downstream scoring in a single end-to-end path. In practice, some runs partially succeed: the library can load, the raw data can be parsed, and a subset of targets can be scored, but the calibration stage may not produce enough significant identifications to safely initialize all downstream components. Existing workflows often degrade poorly in this situation. They either hard-fail the run or continue in a way that mixes calibrated and uncalibrated evidence without a clean confidence model.

The practical gap is not "how to calibrate mass spectrometry data" in general. The gap is how to detect insufficient calibration confidence early, branch the workflow into alternate execution paths, continue with constrained analysis where justified, and emit result objects with explicit confidence partitions that distinguish identification confidence from quantitation confidence and calibration sufficiency.

### Candidate Novelty

The inventive concept is a workflow controller that:

- computes calibration sufficiency metrics during the run rather than treating calibration as only pass/fail
- uses those metrics to route each batch, tranche, precursor family, or peptide cohort into one of multiple scoring modes
- selectively enables or disables downstream steps such as ppm optimization, feature extraction, neural reranking, or protein rollup based on the confidence profile of the upstream calibration state
- emits structured results carrying separate confidence dimensions, for example:
  - identification confidence
  - quantitation confidence
  - calibration confidence
  - fallback-path provenance
- optionally retries subregions of the run with alternate extraction windows, alternate feature sets, or alternative local calibration estimators before declaring a global failure

This is narrower than claiming calibration, recalibration, or optimization themselves. The target claim surface is the orchestration logic, confidence propagation, and structured degraded-mode execution.

### Why This Angle Is Better Than Broad Calibration Claims

Broad DIA analysis, targeted extraction, chromatogram-library style processing, and ML-assisted rescoring already have substantial prior art. The claim should not be "perform DIA search with calibration and scoring." It should be the controller that decides when and how to continue under partial calibration failure, and how to represent the resulting uncertainty in the output.

### Draft Independent Claim

1. A computer-implemented method for analyzing data-independent acquisition mass spectrometry data, the method comprising:
   - receiving a spectral library, a sequence database, configuration parameters, and raw or converted mass spectrometry data;
   - performing an initial calibration procedure on at least a subset of candidate analytes to generate one or more calibration sufficiency metrics;
   - determining, using the calibration sufficiency metrics, that a first subset of candidate analytes satisfies a first confidence condition and that a second subset of candidate analytes satisfies a second, lower confidence condition;
   - routing the first subset through a first downstream analysis path that uses calibration-derived parameters;
   - routing the second subset through a second downstream analysis path that restricts or modifies one or more of mass-accuracy optimization, extraction-window selection, feature generation, discriminant scoring, neural-network reranking, or protein-level aggregation;
   - generating result records that encode separate confidence values for identification and quantitation together with provenance describing which downstream analysis path was applied; and
   - writing the result records to a machine-readable output artifact.

### Dependent Claim Candidates

2. The method of claim 1, wherein the calibration sufficiency metrics include one or more of:
   - number of significant target-decoy separations
   - weighted false-discovery proxy
   - retention-time fit quality
   - mass-error dispersion
   - confidence stability across batches

3. The method of claim 1, wherein the second downstream analysis path disables mass-accuracy optimization and uses conservative extraction windows.

4. The method of claim 1, wherein the second downstream analysis path applies a local calibration model computed for a subset of retention-time bins, m/z bins, ion-mobility bins, or acquisition windows.

5. The method of claim 1, wherein the result record includes a machine-readable field indicating whether a score was produced from a fully calibrated path, a partially calibrated path, or a fallback path.

6. The method of claim 1, wherein protein-level aggregation excludes candidates from the second downstream analysis path unless a quantitation-confidence threshold is satisfied.

7. The method of claim 1, further comprising retrying at least one subset of candidate analytes using an alternate feature set or alternate extraction width before assigning the second downstream analysis path.

### Implementation Hooks in This Repo

- Workflow orchestration: `WorkFlowsLib/src/PythiaDIAFFWorkflow.cpp`
- Calibration path: `WorkFlowsLib/src/PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.cpp`
- Mass-accuracy tuning: `WorkFlowsLib/src/PythiaDIAFFWorkflowAlgos/OptimizeMassAccuracyPPMSettertron.cpp`
- Shared batch processing: `WorkFlowsLib/src/PythiaDIAFFWorkflowAlgos/PythiaDIAFFWorkflowSharedMethods.cpp`

### Claiming Guidance

Keep the independent claim centered on:

- calibration sufficiency estimation
- workflow branching based on that estimate
- different downstream analysis modes
- explicit multi-axis confidence outputs

Avoid claiming:

- generic DIA search
- generic mass recalibration
- generic neural-network reranking

## Draft 2: Interference-Conditioned Local Null Modeling for DIA FDR

### Working Title

Interference-conditioned local null model generation for false-discovery estimation in DIA proteomics.

### Problem Statement

Conventional target-decoy methods and library-decoy generation schemes are useful, but they typically create a global null distribution or a library-level null distribution. DIA data often violate the assumptions behind a single global null because interference structure changes with precursor density, coelution, fragment-sharing, retention-time region, acquisition window, charge state, and sample composition.

The practical gap is that false-discovery estimation is often performed using decoys that match target peptides at a global level, while the actual error process is local and interference-dependent. A local interference-heavy region of the run can behave very differently from a sparse low-conflict region. The same score threshold can therefore be over-conservative in one context and anti-conservative in another.

### Candidate Novelty

The inventive concept is a null-model generator and scorer that:

- measures an interference context for each target candidate or candidate cohort
- constructs a local decoy cohort or null cohort matched to that interference context
- estimates local error statistics using those matched cohorts rather than only a global decoy pool
- optionally blends local and global null estimates according to sample size, confidence, or stability
- outputs context-aware q-values or posterior error estimates that are conditioned on the interference regime of the candidate

The interference context can include:

- fragment-ion sharing with other library entries
- coelution density
- acquisition-window occupancy
- precursor-isotope ambiguity
- local library redundancy
- conflict assays or predicted interference burden

### Why This Angle Is Better Than Broad Decoy Claims

Generic target-decoy analysis and target-like decoy generation are already well developed. The better claim surface is not "generate decoys." It is generating and applying local null cohorts conditioned on measured interference structure, then using those cohorts for candidate-level or cohort-level FDR estimation.

### Draft Independent Claim

1. A computer-implemented method for estimating false-discovery rates in data-independent acquisition proteomics, the method comprising:
   - receiving a set of target candidates derived from a spectral library and mass spectrometry data;
   - determining, for at least one target candidate, an interference context comprising one or more measurements of fragment conflict, precursor density, coelution burden, acquisition-window occupancy, or library redundancy;
   - generating or selecting a matched null cohort for the target candidate based on the interference context;
   - computing a context-conditioned error estimate for the target candidate using scores from the matched null cohort; and
   - reporting an identification confidence metric derived from the context-conditioned error estimate.

### Dependent Claim Candidates

2. The method of claim 1, wherein the matched null cohort is generated by selecting decoy candidates that match the target candidate with respect to:
   - precursor m/z range
   - charge state
   - retention-time neighborhood
   - fragment-sharing burden
   - acquisition-window membership

3. The method of claim 1, wherein the interference context further includes a metric of predicted chromatographic overlap or observed fragment-ion competition.

4. The method of claim 1, wherein the context-conditioned error estimate is blended with a global target-decoy estimate according to a confidence weight determined from cohort size or score stability.

5. The method of claim 1, wherein different target candidates within a single run are assigned different null cohorts based on different interference contexts.

6. The method of claim 1, wherein the reported identification confidence metric is accompanied by a separate quantitation-reliability metric.

7. The method of claim 1, wherein the matched null cohort is stored or serialised as provenance linked to the corresponding target candidate.

### Implementation Hooks in This Repo

- Target/decoy construction: `AlgorithmsFFLib/src/TargetDecoyCandidatePairManager.h`
- Candidate scoring: `AlgorithmsFFLib/src/TargetDecoyCandidatePairScoretron.h`
- Discriminant scoring: `AlgorithmsFFLib/src/DiscriminantScoretron.h`
- Neural reranking: `AlgorithmsFFLib/src/FDRCLassifierNeuralNet.h`
- Output serialization: `FileReadersLib/src/ParquetReader.h`

### Claiming Guidance

Keep the independent claim centered on:

- measuring interference context
- selecting or generating a matched local null cohort
- computing context-conditioned FDR or related confidence values

Avoid claiming:

- generic target-decoy FDR
- generic library decoys
- generic machine-learning rescoring without the local-null conditioning step

## Suggested Next Step for Counsel

If these move forward, the next useful deliverable is not code. It is a prior-art matrix with columns for:

- claim element
- likely academic prior art
- likely patent prior art
- novelty argument
- fallback narrowing language

That matrix will show whether the inventive center is actually in the controller logic, the local null construction, or the output confidence model.
