# Radiant DIA API Reference

This document describes the public-facing interfaces exposed by the `radiant` repository. It focuses on the command-line surface, configuration schema, and the main C++ classes intended to be composed by the shipped applications.

## Scope

This reference covers:

- CLI entry points
- config-file schema through `PythiaParameters`
- core workflow and reader classes
- scoring and reranking classes
- output contracts

It does **not** attempt to document every internal helper in `AlgorithmsFFLib`. For internals such as `CandidateScorertron`, `TurboXIC`, or `XICPeakManager`, this guide treats them as implementation details unless they are surfaced through public workflow APIs.

## Conventions

### Error model

Most library methods return `Error::Err`:

```cpp
enum Err {
    eNoError = 0,
    eError,
    eEmptyContainerError,
    eValueError,
    eMemoryError,
    eFileError,
    eFileIncorrectTypeError,
    eNetworkError,
    eFunctionNotImplemented
};
```

Implications:

- `eNoError` means success
- non-zero values are propagated aggressively via macros such as `ree`
- methods generally log the file and line on failure before returning

### Type system

The public APIs rely heavily on Qt containers and strings:

- `QString`
- `QStringList`
- `QVector<T>`
- `QMap<K, V>`
- `QHash<K, V>`
- `QSharedPointer<T>`

### Ownership model

Ownership varies by class:

- workflow objects such as `PythiaDIAFFWorkflow` own internal mutable state and should be treated as single-run objects
- `MsReaderPointerAcc` exposes the active reader as `QSharedPointer<MsReaderBase> ptr`
- `FDRCLassifierNeuralNet` allocates `CandidateClassifier` instances internally

## CLI API

## `RadiantDIA`

Primary executable for DIA library search and reranking.

Usage:

```text
RadiantDIA [options] fraglibff-path fasta-path config-path data-file-path
```

Arguments:

- `fraglibff-path`: fragment library file
- `fasta-path`: FASTA file
- `config-path`: `.radiantConfig`, `.pythiaConfig`, or `.toml`
- `data-file-path`: `.prqFF`, `.mzML`, or Bruker `.d`

Options:

- `--output-folder <path>`: write the `.radiantDIA` result into a specific folder

Primary code:

- `Apps/RadiantDIA/main.cpp`
- `Apps/RadiantDIA/src/CommandLineParser.cpp`

## `LycaonMzMLTransmorgitron`

Auxiliary conversion executable.

Usage:

```text
LycaonMzMLTransmorgitron mzML-path
```

Purpose:

- accepts a single `.mzML` file
- fronts the mzML-to-Parquet conversion workflow

Primary code:

- `Apps/LycaonMzMLTransmorgitron/main.cpp`
- `Apps/LycaonMzMLTransmorgitron/src/CommandLineParser.cpp`

## `FeatureOptimizerPrime`

Auxiliary optimization executable.

Usage shape:

```text
FeatureOptimizerPrime fraglibff-path fasta-path pythia-path data-file-path
```

Purpose:

- runs optimization-oriented workflow paths against the same core data contracts used by `RadiantDIA`

Primary code:

- `Apps/FeatureOptimizerPrime/main.cpp`
- `Apps/FeatureOptimizerPrime/src/CommandLineParser.cpp`

## Configuration API

Configuration is represented in code by `PythiaParameters` and parsed by `PythiaParameterReader`.

### Parser entry points

```cpp
static PythiaParameters genericPythiaParametersForTests();
static Err buildPythiaParameters(
    const QString &pythiaParametersFilePath,
    PythiaParameters *pythiaParameters
);
```

### Supported config file extensions

- `.radiantConfig`
- legacy `.pythiaConfig`
- `.toml`

### Config sections and keys

## `[General]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `threadCount` | `int` | system CPU count if unset/non-positive | worker count |
| `verbosity` | `int` | `0` | log verbosity |
| `writeRadiantDIA` | `bool` | `true` | write final `.radiantDIA` output |
| `useLazyLoading` | `bool` | `false` | enable on-demand mz-target loading where supported |
| `reannotate` | `bool` | `false` | re-map peptides to protein groups from FASTA |
| `shortReport` | `bool` | `false` | write truncated result rows |

## `[LibraryParams]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `chargeStateMin` | `int` | required via validity | minimum precursor charge |
| `chargeStateMax` | `int` | required via validity | maximum precursor charge |
| `mzMinMS2` | `double` | `300.0` | lower fragment-library m/z bound |
| `mzMaxMS2` | `double` | `1999.0` | upper fragment-library m/z bound |
| `peptideLengthMin` | `int` | `7` | minimum peptide length |
| `peptideLengthMax` | `int` | `30` | maximum peptide length |
| `trancheSizeMax` | `int` | `20000` | tranche size for parallel scoring |
| `useAlternativeDecoys` | `bool` | `false` | alternate decoy fragment-shift behavior |

## `[MS1Params]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `precursorExtractionWindowThomsons` | `double` | `0.0` | precursor isolation margin |
| `ms1ExtractionWidthPPM` | `double` | `20.0` | baseline MS1 extraction width |
| `ms1ExtractionWidthPPMOverride` | `double` | `-1.0` | force MS1 extraction width, bypassing learned value |

## `[MS2Params]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `filterLengthIntegration` | `int` | `5` | integration filter length |
| `filterLengthMS2` | `int` | `3` | MS2 smoothing/filter length |
| `ionsSharedToReject` | `int` | `4` | interference rejection threshold |
| `ms2ExtractionWidthPPM` | `double` | `20.0` | baseline MS2 extraction width |
| `ms2ExtractionWidthPPMOverride` | `double` | `-1.0` | force MS2 extraction width |
| `minMs2FragCount` | `int` | `2` | minimum fragment support |
| `scanTimeWindowStDevs` | `float` | `3` | RT extraction window in standard deviations |
| `subtractShadows` | `bool` | `true` | remove interfering shadow signal |
| `smoothCountMS2` | `int` | `1` | MS2 smoothing passes |
| `stopThresholdFractionMS2` | `float` | `0.65` | integration cutoff fraction |
| `peakDifferenceThresholdFraction` | `float` | `0.1` | peak-difference guard |
| `calibrationTrainingVolume` | `int` | `1000` | calibration scoring volume |
| `peakCenter` | `int` | `-1` | optional forced peak center |
| `topNIntegrations` | `int` | `10` | number of fragment integrations retained |
| `maxAnchorColumnIndex` | `int` | `6` | anchor fragment index ceiling |
| `rtBinning` | `int` | `20` | RT binning for some optimization paths |

## `[FdrParams]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `percentFDR` | `double` | `1.0` | target FDR percentage |
| `reportDecoys` | `bool` | `false` | keep/report decoys in final output |

## `[PeakIntegrationParams]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `filterLength` | `int` | unset until config | integration filter width |
| `smoothCount` | `int` | unset until config | integration smoothing count |
| `sigma` | `double` | unset until config | smoothing sigma |
| `signalToNoiseRatio` | `double` | `2.0` | signal/noise threshold |
| `stopThresholdFraction` | `float` | `0.65` | stop threshold for generic integration |

## `[FeatureFinderParams]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `minScanCount` | `int` | `3` | minimum hill/feature scan span |
| `skipScanCount` | `int` | `2` | scan skip interval for feature finding |

## `[NeuralNetParams]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `epochs` | `int` | `12` | NN training epochs |
| `baggingSize` | `int` | `12` | ensemble size; clamped to `>= 2` |
| `learningRate` | `float` | `0.003` | optimizer learning rate |
| `nodesFraction` | `double` | `0.5` | hidden-layer width multiplier |
| `focalLossGamma` | `float` | `0.0` | focal-loss gamma; non-positive behaves like BCE |
| `parallelNeuralNets` | `bool` | `false` | intended NN training parallelism toggle |

### Config caveat

The parser reads `writeRadiantDIA`, not `writePythiaDIA`. One bundled sample config still uses the legacy `writePythiaDIA` key, so that key currently falls back to the default rather than being explicitly parsed.

## Core Workflow API

## `ConvertMzMLToParquetWorkFlow`

Header:

```cpp
class ConvertMzMLToParquetWorkFlow {
public:
    static Err convertMzMLToParquet(
        const QString &mzmlFilePath,
        QString *outputFilePath
    );
};
```

Purpose:

- convert mzML input into the Parquet-based `.prqFF` representation used by faster local paths

Use when:

- precomputing a Parquet representation before full search
- normalizing inputs for repeated local testing

## `PythiaDIAFFWorkflow`

Header:

```cpp
Err init(
    const PythiaParameters &pythiaParameters,
    const QString &fragLibUri,
    const QString &fastaUri,
    const QString &outputFolderPath
);

Err processFile(const QString &msDataFilePath);

Err setCalibratomaticFeatures(const QVector<Features> &features);
Err setPPMOptimizationFeatures(const QVector<Features> &features);
Err setNeuralNetFeatures(const QVector<Features> &features);

[[nodiscard]] ResultsSummary resultsSummary() const;
```

Responsibilities:

- own one end-to-end Radiant DIA search run
- sequence all stages from library load through final output
- manage mutable state for:
  - target/decoy candidate pairs
  - calibration state
  - learned weights
  - candidate score collections
  - result summaries

Important behavior:

- `init()` loads library rows and builds candidate pairs
- `processFile()` opens the MS data file, calibrates, optimizes, scores, reranks, annotates, and writes results
- `set*Features()` lets callers replace the default calibration, optimization, and neural-net feature sets before running

## `TargetDecoyCandidatePairManager`

Header:

```cpp
Err init(
    const PythiaParameters &pythiaParameters,
    QList<FragLibReaderRow> *fragLibReaderRows
);

bool isInit() const;
int targetsCount() const;

Err getTargetDecoyCandidatePairPointers(
    QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairsPntrs
);

static Err peptideStringWithModsFromPeptideSequenceChargeKey(
    const PeptideSequenceChargeKey &peptideSequenceChargeKey,
    PeptideStringWithMods *peptideStringWithMods,
    int *charge
);
```

Responsibilities:

- turn fragment-library rows into target/decoy search candidates
- enforce peptide length and charge constraints
- apply decoy-generation rules
- deduplicate and flag problematic target/decoy collisions

Use when:

- you need the search space without running the full workflow

## `TargetDecoyCandidatePairScoretron2`

Header:

```cpp
Err init(
    const PythiaParameters &pythiaParameters,
    MsReaderPointerAcc *msReaderPointerAcc
);

Err scoreTargetDecoyPairs(
    const QVector<Features> &features,
    int topNMS2Ions,
    const MsCalibratomatic &msCalibratomatic,
    float minPeakCount,
    int threadCount,
    bool useTopNIntegrationsParameter,
    const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXicPntrs,
    const QVector<float> &weights,
    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> *candidateScoresPairsVec
);

Err scoreTargetDecoyPairs(
    const QVector<Features> &features,
    int topNMS2Ions,
    const MsCalibratomatic &msCalibratomatic,
    float minPeakCount,
    int threadCount,
    bool useTopNIntegrationsParameter,
    const QVector<MsScanInfo> &msScanInfos,
    const QVector<float> &weights,
    QVector<TargetDecoyCandidatePair*> *targetDecoyCandidateAllPntrs,
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> *candidateScoresPairsVec
);

bool isInit() const;
Err setPythiaParameters(const PythiaParameters &pythiaParameters);
QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>>* diaTargetFrames();
QMap<ScanNumber, ScanPoints>* ms1ScanNumberVsScanPoints();
QMap<MzTargetKey, MsFrame*> mzTargetKeyVsMsFramePntr();
Err reloadTurboXICMS1();
Err buildMzTargetKeyVsMsFrames();
```

Responsibilities:

- bridge from loaded MS data to scored target/decoy candidates
- own MS1 frame state, per-window MS2 frame state, and averagine tables
- prepare parallel scoring inputs
- execute the heavy feature-extraction and scoring phase

Use when:

- building custom search flows outside the stock workflow
- rescoring candidate sets with different feature selections or tolerances

## Reader and Format API

## `MsReaderPointerAcc`

Header:

```cpp
Err openFile(const QString &filePath);
Err openFile(
    const QString &filePath,
    const QString &columnToFilterBy,
    const QPair<double, double> &filterRange
);
Err openFile(
    const QString &filePath,
    const QString &columnToFilterBy
);

bool isInit();
void setUseLazyLoading(bool useLazyLoading);
bool useLazyLoading() const;

QSharedPointer<MsReaderBase> ptr;
```

Responsibilities:

- choose and initialize the correct concrete MS reader for the input path
- expose the active reader through the shared `ptr`
- support lazy-loading mode where reader implementations provide it

Supported file families:

- `.mzML`
- `.prqFF`
- `.d`

## `FragLibReader`

Header:

```cpp
static Err getFragLibReaderRows(
    const QString &fragLibFilePath,
    bool useAlternativeDecoys,
    QList<FragLibReaderRow> *fragLibReaderRows
);
```

Responsibilities:

- normalize supported fragment-library formats into `FragLibReaderRow` objects

## `PythiaParameterReader`

Header:

```cpp
static PythiaParameters genericPythiaParametersForTests();
static Err buildPythiaParameters(
    const QString &pythiaParametersFilePath,
    PythiaParameters *pythiaParameters
);
```

Responsibilities:

- parse config files into the in-memory parameter object used by all major workflow stages

## `ParquetReaderInputBase`

Key extension points:

```cpp
virtual QMap<QString, QVariant> map();
virtual Err initFromRead(const ParquetReaderInputBase &row);

template<typename T>
static QVector<QSharedPointer<ParquetReaderInputBase>> convertInputStructToSharedPointers(...);

template<typename T>
static Err convertSharedPointersToInputStruct(...);
```

Purpose:

- define the row contract used by Parquet-backed result and intermediate artifacts
- support serialization/deserialization of derived row types such as `CandidateScoresReaderRow`

## `ParquetReader`

Purpose:

- generic Parquet read/write adapter for row types derived from `ParquetReaderInputBase`

Use when:

- writing result artifacts
- reading/writing intermediate tables or cached representations

## Scoring and Model API

## `DiscriminantScoretron`

Header:

```cpp
static QVector<Features> featuresCalibration();
static QVector<Features> featuresOptimization();
static QVector<Features> featuresNeuralNetwork();

static Err trainLDAClassifier(
    const QVector<QPair<FeaturesArrayDecoys*, FeaturesArrayDecoys*>> &targetDecoyCandidateScoresPair,
    int verbosity,
    QVector<float> *weights
);

static Err applyWeights(
    const QVector<float> &weights,
    int threadCount,
    const QVector<FeaturesArray*> &featureArrayPntrs,
    QVector<float> *discriminantScores
);

static QVector<float> defaultWeights(const QVector<Features> &features);

static Err convertScoreCandidatesToFeaturesArrays(
    const QVector<Features> &features,
    const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScoresTargetVsDecoy,
    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy
);
```

Responsibilities:

- define the feature subsets used in each workflow phase
- train and apply the linear discriminant stage
- convert candidate scores into feature arrays suitable for classification

## `FDRCLassifierNeuralNet`

Header:

```cpp
Err init(
    int epochs,
    int baggingSize,
    int batchSize,
    double learningRate,
    double nodesFraction,
    float focalLossGamma,
    int threadCount
);

Err exec(
    const QVector<QVector<float>> &xData,
    const QVector<float> &yData,
    int seed,
    int verbosity,
    QVector<float> *meanPredictions
);

Err trainClassifier(
    const QVector<QVector<float>> &xData,
    const QVector<float> &yData,
    int seed,
    int verbosity
);

Err predictBaggedClassifiers(
    const QVector<QVector<float>> &allDataVecs,
    QVector<float> *meanPredictions
);

static Err countScoreCandidatesByFDR(...);
static Err outputFDRResults(...);
static Err outPutFDRCounts(...);
```

Responsibilities:

- train the bagged neural-network ensemble
- average predictions across the ensemble
- compute and summarize q-value/FDR counts

Important behavior:

- bagging size is clamped to at least `2`
- `exec()` trains and predicts in one call

## `CandidateClassifier`

Header:

```cpp
bool trainCandidateClassifier(
    const QVector<QVector<float>> &xData,
    const QVector<float> &yData,
    int epochsMax,
    int batchSize,
    double learningRate,
    int seed,
    double nodeFraction,
    float focalLossGamma,
    int verbosity
) const;

bool predict(
    const QVector<QVector<float>> &xData,
    QVector<float> *predictions
) const;
```

Responsibilities:

- encapsulate the actual libtorch training and inference logic used by `FDRCLassifierNeuralNet`

Important distinction:

- this class returns `bool`, not `Err`
- it is a lower-level ML adapter rather than a full workflow component

## Results and Output Contracts

## Primary result file

Radiant's main output file has extension:

- `.radiantDIA`

It is written as Parquet and contains either:

- `CandidateScoresReaderRow`
- `CandidateScoresReaderRowTrunc`

depending on `shortReport`.

## `ResultsSummary`

The workflow exposes a lightweight summary object:

```cpp
struct ResultsSummary {
    int altPSMCountNeuralNet;
    int psmCountNeuralNet;
    int decoyCountNeuralNet;
    int entrapCountNeuralNet;
    int psmCountCorrectedFMRNeuralNet;
    float eFDRNeuralNet;

    int psmCountLDA;
    int decoyCountLDA;
    int entrapCountLDA;
    int psmCountCorrectedFMRLDA;
    float eFDRLDA;
};
```

This is populated mainly when the protein reannotation and entrapment-report logic runs.

## End-to-End Integration Sequence

The intended calling pattern for the main search surface is:

```cpp
PythiaParameters params;
Err e = PythiaParameterReader::buildPythiaParameters(configPath, &params);

PythiaDIAFFWorkflow workflow;
e = workflow.init(params, fragLibPath, fastaPath, outputFolderPath);
e = workflow.processFile(msDataFilePath);

ResultsSummary summary = workflow.resultsSummary();
```

For conversion-only paths:

```cpp
QString outputFilePath;
Err e = ConvertMzMLToParquetWorkFlow::convertMzMLToParquet(mzmlPath, &outputFilePath);
```

## Known API Caveats

## Legacy naming

Many public types still use legacy `Pythia` naming:

- `PythiaParameters`
- `PythiaDIAFFWorkflow`
- helper methods named `writePythiaDIA`

Those are current code symbols, but product-facing docs should map them mentally to `Radiant DIA`.

## Sample config drift

The parser expects `writeRadiantDIA`; one bundled example still uses `writePythiaDIA`. That sample key is stale.

## Mutable workflow state

`PythiaDIAFFWorkflow` accumulates mutable state across initialization and scoring. Treat it as a per-run object rather than a stateless service.

## End-to-end fixture reliability

The bundled test assets are enough to boot and exercise the binary, but they are not guaranteed to produce a successful `.radiantDIA` output under the current defaults. Calibration and downstream scoring can still fail to produce a complete run.
