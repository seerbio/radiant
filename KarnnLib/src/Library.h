/*
Copyright 2020, Vadim Demichev

This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license,
visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

Uncomment "#define _NO_THERMORAW" in RAWReader.cpp, MSReader.cpp and MSReader.h to enable building without Thermo MSFileReader installed
*/

#ifndef PYTHIADIACPP_TEST
#define PYTHIADIACPP_TEST

#include "KarnnLib_Exports.h"

#include "cpu_info.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <string>
#include <regex>
#include <algorithm>
#include <random>
#include <atomic>
#include <thread>
#include <math.h>
#include <list>
#include <set>
#include <map>

#include <QDebug>

#define Min(x, y) ((x) < (y) ? (x) : (y))
#define Max(x, y) ((x) > (y) ? (x) : (y))
#define Abs(x) ((x) >= 0 ? (x) : (-(x)))
#define Sgn(x) ((x) >= 0.0 ? (1) : (-1))
#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * Sqr(x))

#define Throw(x) { qDebug() << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << x << "\n"; exit(-1); }
#define Warning(x) { if (Verbose >= 1) qDebug() << "WARNING: " << x << "\n"; }

typedef std::chrono::high_resolution_clock Clock;
static auto StartTime = Clock::now();
#define Minutes(x) (floor((x) * 0.001 * (1.0/60.0)))
#define Seconds(x) (floor((x) * 0.001) - 60 * Minutes(x))
inline void Time() {
    double time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - StartTime).count();
    int min = Minutes(time);
    int sec = Seconds(time);
    qDebug() << "[" << min << ":";
    if (sec < 10) qDebug() << "0";
    qDebug() << sec << "] ";
}

const double E = 0.000000001;
const double INF = 10000000.0;

int Threads = 1;
int iN = 12;
int nnIter = 11;

bool LCAllScores = false;

bool Standardise = true;
float StandardisationScale = 1.0;

int QuantMode = 0; // 0 - full peak integration, 1 - fixed window
double PeakBoundary = 5.0;
bool NoIfsRemoval = false;
bool NoFragmentSelectionForQuant = false;
bool QuantFitProfiles = false;
bool QuantAdjacentBins = true;
bool RestrictFragments = false;

bool BatchMode = true;
int MinBatch = 2000;
int MaxBatches = 10000;
int Batches = 1;

const int MaxLibSize = 1000000;

bool RTProfiling = false;

int ScanFactor = 150;
int WindowRadius = 0;
double ScanScale = 2.2;
bool InferWindow = true;
bool IndividualWindows = false;
bool Calibrate = true;

const int MaxCycleLength = 2000;
double MinPeakHeight = 0.01;

int iRTTargetTopPrecursors = 50;
int iRTTopPrecursors;
int iRTRefTopPrecursors = 250;
double iRTMaxQvalue = 0.1; // only precursors with lower q-value are used for iRT estimation
double MassCalQvalue = 0.1;
double MassCalMs1Corr = 0.90;
int MassCalBinsMax = 1;
int MassCalBins = 1;
int MassCalBinsMs1 = 1;
int MinMassCalBin = 500;
bool MassCalFilter = true;

int RTSegments = 20;
int MinRTPredBin = 20;

bool ForceScanningSWATH = false;
bool ForceNormalSWATH = false;
bool UseQ1 = true;
bool ForceQ1 = false;
bool Q1Cal = false;
bool Q1CalLinear = false;
int MaxQ1Bins = 3;
int MaxIfsSpan = 4;
bool ForceQ1Apex = false;

bool RefCal = false;
int RTWindowIter = 4;
const int RTRefWindowFactor = 10;
const int RTWindowFactor = 40;
const double RTWindowLoss = 0.20;
const double RTWindowMargin = 2.0;
double MinRTWinFactor = 1.0;
int RTWinSearchMinRef = 10;
int RTWinSearchMinCal = 100;

int nnBagging = 12;
int nnEpochs = 1;
int nnHidden = 5;
double nnLearning = 0.003;
double Regularisation = E;
bool nnFilter = true;
bool nnCrossVal = false;

bool TightMassAcc = true;
bool TightMassAauxForCal = false;
double TightMassAccRatioOne = 0.45;
double TightMassAccRatioTwo = 0.2;

bool UseIsotopes = true;
int IDsInterference = 1;
double InterferenceCorrMargin = 3.0;
bool StrictIntRemoval = false;
bool MS2Range = true;
bool GuideClassifier = false;

bool RTWindowedSearch = true;
bool DisableRT = false;

const int CalibrationIter = 4;
bool ForceMassAcc = false;
bool IndividualMassAcc = false;
double CalibrationMassAccuracy = 100.0 / 1000000.0;
double GlobalMassAccuracy = 20.0 / 1000000.0;
double GlobalMassAccuracyMs1 = 20.0 / 1000000.0;
double GeneratorAccuracy = 5.0 / 1000000.0;
double SptxtMassAccuracy = 5.0 / 1000000.0;
double MaxProfilingQvalue = 0.001;
double MaxQuantQvalue = 0.01;
double ProteinQuantQvalue = 0.01;
double ProteinIDQvalue = 0.01;
double ProteinQCalcQvalue = 0.02; // must be 0.01 or greater
double PeakApexEvidence = 0.99; // force peak apexes to be close to the detection evidence maxima
double MinCorrScore = 0.5;
double MinMs1Corr = -INF;
double MaxCorrDiff = 2.0;
bool ForceMs1 = false;
bool MS1PeakSelection = true;
bool MaxLFQ = true;
int TopN = 1;

bool TranslatePeaks = false; // translate peaks between precursors belonging to the same elution group
bool ExcludeSharedFragments = true; // exclude from quantification fragments that are shared between light and heavy labelled peptides

int MinMassDeltaCal = 20;
int MinMassDeltaCalRef = 8;

int MinCalRec = 100;
int MinCal = 1000;
int MinClassifier = 2000;

double FilterFactor = 1.5;

bool Convert = false;
bool MGF = false;
bool Reannotate = false;
bool UseRTInfo = false; // use .quant files created previously for RT profiling
bool UseQuant = false; // use .quant files created previously; implies UseRTInfo
bool QuantOnly = false; // quantification will be performed anew using identification info from .quant files
bool ReportOnly = false; // generate report from quant files
bool GenRef = false; // update the .ref file with newly computed data
std::string args, lib_file, learn_lib_file, out_file = "report.tsv", prosit_file = "lib.prosit.csv", predictor_file = "lib.predicted.speclib", out_lib_file = "lib.tsv", out_dir = "", out_gene_file = "report.genes.tsv", temp_folder = "";
std::string ref_file, gen_ref_file, all_fastas;
std::vector<std::string> ms_files, fasta_files, fasta_filter_files;
std::set<std::string> failed_files;
bool FastaSearch = false;
bool GenSpecLib = false;
bool ProfileQValue = true;
bool RTLearnLib = false;
bool SaveOriginalLib = false;
bool SaveCalInfo = false;
bool iRTOutputLibrary = true;
bool PredictorSaved = false;
bool LibrarySaveRT = false;

int Verbose = 1;
bool ExportWindows = false;
bool ExportLibrary = false;
bool ExportDecoys = false;
bool GuideLibrary = false;
bool OverwriteLibraryPGs = true;
bool InSilicoRTPrediction = true;
bool ReverseDecoys = false;
bool ForceFragRec = false;
bool ExportRecFragOnly = true;
bool GenFrExclusionInfo = false;

bool ExportProsit = false;

int VisWindowRadius = 8;
std::vector<std::string> Visualise;
bool SaveXICs = false;
int XICRadius = 20;

enum {
    loss_none, loss_H2O, loss_NH3, loss_CO,
    loss_N, loss_other
};

std::string getLoss(int lossIndex) {

    switch (lossIndex) {

        case 0:
            return "";
        case 1:
            return "-H2O";
        case 2:
            return "-NH3";
        case 3:
            return "a";
        case 4:
            return "N";
        case 5:
            return "other";
        default:
            return "";
    }

    return "";
}


const double Loss[loss_other + 1] = { 0.0, 18.011113035, 17.026548, 27.994915, 0.0, 0.0 };

int MaxRecCharge = 19;
int MaxRecLoss = loss_N;

double MinMs1RangeOverlap = 0.25;
double QuadrupoleError = 1.0;

std::vector<int16_t> Cut = { ('K' << 8) | '*', ('R' << 8) | '*' }, NoCut;

int MissedCleavages = 1;
int MinPeptideLength = 7;
int MaxPeptideLength = 30;
bool NMetExcision = false;
bool AddLosses = false;
bool OriginalModNames = false;
int MaxVarMods = 1;
int MinFrAAs = 3;
double MinFrMz = 200.0, MaxFrMz = 1800.0;
double MinPrMz = 300.0, MaxPrMz = 1800.0;
double MinGenFrInt = 0.00001;
double MinRelFrHeight = 0.00001;
double MinGenFrCorr = 0.5;
double MinRareFrCorr = 0.95;
double FrCorrBSeriesPenalty = 0.1;
int MinPrCharge = 1;
int MaxPrCharge = 10;
int MinGenFrNum = 4;
int MinOutFrNum = 4;
int MinSearchFrNum = 3;
const int MaxGenCharge = 4;

double QFilter = 1.0;
double RubbishFilter = 0.5;
int MinRubbishFactor = 5;
int MinNonRubbish = 5000;

bool InferPGs = true;
bool PGsInferred = false;

double ReportQValue = 0.01;
double ReportProteinQValue = 1.0;
double MatrixQValue = 0.01;
int PGLevel = 2, PGLevelSet = PGLevel; // 0 - ids, 1 - names, 2 - genes
std::vector<std::string> ImplicitProteinGrouping = { "isoform IDs", "protein names", "genes" };
bool IndividualReports = false;
bool SwissProtPriority = true;
bool ForceSwissProt = false;
bool FastaProtDuplicates = false;
bool SpeciesGenes = false;

bool ExtendedReport = true;
bool ReportLibraryInfo = false;
bool DecoyReport = false;
bool SaveCandidatePSMs = false;
bool MobilityReport = false;
bool SaveMatrices = false;

bool RemoveQuant = false;
bool QuantInMem = false;

#define Q1 true
#define AAS true
#define LOGSC false
#define REPORT_SCORES false
#define EXTERNAL 1
#define ELUTION_PROFILE false
const int ElutionProfileRadius = 10;

enum {
    libPr, libCharge, libPrMz,
    libiRT, libFrMz, libFrI,
    libPID, libPN, libGenes, libPT, libIsDecoy,
    libFrCharge, libFrType, libFrNumber, libFrLoss,
    libQ, libEG, libFrExc, libCols
};

std::vector<std::string> library_headers = {
        " ModifiedPeptide LabeledSequence FullUniModPeptideName \"FullUniModPeptideName\" modification_sequence ",
        " PrecursorCharge Charge \"PrecursorCharge\" prec_z ",
        " PrecursorMz \"PrecursorMz\" Q1 ",
        " iRT iRT RetentionTime NormalizedRetentionTime Tr_recalibrated \"Tr_recalibrated\" RT_detected ",
        " FragmentMz ProductMz \"ProductMz\" Q3 ",
        " RelativeIntensity RelativeFragmentIntensity RelativeFragmentIonIntensity LibraryIntensity \"LibraryIntensity\" relative_intensity ",
        " UniprotId UniProtIds UniprotID \"UniprotID\" uniprot_id ",
        " Protein Name Protein.Name ProteinName \"ProteinName\" ",
        " Genes Gene Genes GeneName \"Genes\" ",
        " IsProteotypic Is.Proteotypic Proteotypic ",
        " Decoy decoy \"Decoy\" \"decoy\" ",
        " FragmentCharge FragmentIonCharge ProductCharge ProductIonCharge frg_z \"FragmentCharge\" ",
        " FragmentType FragmentIonType ProductType ProductIonType frg_type \"FragmentType\" ",
        " FragmentNumber frg_nr FragmentSeriesNumber \"FragmentSeriesNumber\" ",
        " FragmentLossType FragmentIonLossType ProductLossType ProductIonLossType \"FragmentLossType\" ",
        " Q.Value QValue Qvalue qvalue \"Q.Value\" \"QValue\" \"Qvalue\" \"qvalue\" ",
        " ElutionGroup ModifiedSequence ",
        " ExcludeFromAssay ExcludeFromQuantification "
};

const double proton = 1.007825035;
const double OH = 17.003288;
const double C13delta = 1.003355;

std::vector<float> Extract;

std::vector<std::pair<std::string, std::string> > UniMod;
std::vector<std::pair<std::string, int> > UniModIndex;
std::vector<int> UniModIndices;

struct MOD {
    std::string name;
    std::string aas;
    float mass = 0.0;
    int label = 0;

    MOD() {};
    MOD(const std::string &_name, float _mass, int _label = 0) { name = _name, mass = _mass; label = _label; }
    void init(const std::string &_name, const std::string &_aas, float _mass, int _label = 0) { name = _name, aas = _aas, mass = _mass; label = _label;  }
    friend bool operator < (const MOD &left, const MOD &right) { return left.name < right.name || (left.name == right.name && left.mass < right.mass); }
};
std::vector<MOD> FixedMods, VarMods;

std::vector<MOD> Modifications = {
        MOD("UniMod:4", (float)57.021464),
        MOD("Carbamidomethyl (C)", (float)57.021464),
        MOD("Carbamidomethyl", (float)57.021464),
        MOD("CAM", (float)57.021464),
        MOD("+57", (float)57.021464),
        MOD("+57.0", (float)57.021464),
        MOD("UniMod:26", (float)39.994915),
        MOD("PCm", (float)39.994915),
        MOD("UniMod:5", (float)43.005814),
        MOD("Carbamylation (KR)", (float)43.005814),
        MOD("+43", (float)43.005814),
        MOD("+43.0", (float)43.005814),
        MOD("CRM", (float)43.005814),
        MOD("UniMod:7", (float)0.984016),
        MOD("Deamidation (NQ)", (float)0.984016),
        MOD("Deamidation", (float)0.984016),
        MOD("Dea", (float)0.984016),
        MOD("+1", (float)0.984016),
        MOD("+1.0", (float)0.984016),
        MOD("UniMod:35", (float)15.994915),
        MOD("Oxidation (M)", (float)15.994915),
        MOD("Oxidation", (float)15.994915),
        MOD("Oxi", (float)15.994915),
        MOD("+16", (float)15.994915),
        MOD("+16.0", (float)15.994915),
        MOD("Oxi", (float)15.994915),
        MOD("UniMod:1", (float)42.010565),
        MOD("Acetyl (Protein N-term)", (float)42.010565),
        MOD("+42", (float)42.010565),
        MOD("+42.0", (float)42.010565),
        MOD("UniMod:255", (float)28.0313),
        MOD("AAR", (float)28.0313),
        MOD("UniMod:254", (float)26.01565),
        MOD("AAS", (float)26.01565),
        MOD("UniMod:122", (float)27.994915),
        MOD("Frm", (float)27.994915),
        MOD("UniMod:1301", (float)128.094963),
        MOD("+1K", (float)128.094963),
        MOD("UniMod:1288", (float)156.101111),
        MOD("+1R", (float)156.101111),
        MOD("UniMod:27", (float)-18.010565),
        MOD("PGE", (float)-18.010565),
        MOD("UniMod:28", (float)-17.026549),
        MOD("PGQ", (float)-17.026549),
        MOD("UniMod:526", (float)-48.003371),
        MOD("DTM", (float)-48.003371),
        MOD("UniMod:325", (float)31.989829),
        MOD("2Ox", (float)31.989829),
        MOD("UniMod:342", (float)15.010899),
        MOD("Amn", (float)15.010899),
        MOD("UniMod:1290", (float)114.042927),
        MOD("2CM", (float)114.042927),
        MOD("UniMod:359", (float)13.979265),
        MOD("PGP", (float)13.979265),
        MOD("UniMod:30", (float)21.981943),
        MOD("NaX", (float)21.981943),
        MOD("UniMod:401", (float)-2.015650),
        MOD("-2H", (float)-2.015650),
        MOD("UniMod:528", (float)14.999666),
        MOD("MDe", (float)14.999666),
        MOD("UniMod:385", (float)-17.026549),
        MOD("dAm", (float)-17.026549),
        MOD("UniMod:23", (float)-18.010565),
        MOD("Dhy", (float)-18.010565),
        MOD("UniMod:129", (float)125.896648),
        MOD("Iod", (float)125.896648),
        MOD("Phosphorylation (ST)", (float)79.966331),
        MOD("UniMod:21", (float)79.966331),
        MOD("+80", (float)79.966331),
        MOD("+80.0", (float)79.966331),
        MOD("UniMod:259", (float)8.014199, 1),
        MOD("Lys8", (float)8.014199, 1),
        MOD("UniMod:267", (float)10.008269, 1),
        MOD("Arg10", (float)10.008269, 1),
        MOD("UniMod:268", (float)6.013809, 1),
        MOD("UniMod:269", (float)10.027228, 1)
};

inline char to_lower(char c) { return c + 32; }

double to_double(const char * s) { // precision ~17 sig figures, locale-independent
    const long long limit = (long long)1 << 58;
    long long tot = 0, sep;
    int exponent = 0;
    bool neg = false;
    const double exps[35] = {
            10.0, 100.0, 1000.0, 10000.0, 100000.0,
            1000000.0, 10000000.0, 100000000.0, 1000000000.0, 10000000000.0,
            100000000000.0, 1000000000000.0, 10000000000000.0, 100000000000000.0, 1000000000000000.0,
            10000000000000000.0, 100000000000000000.0, 1000000000000000000.0, 10000000000000000000.0, 100000000000000000000.0,
            1000000000000000000000.0, 10000000000000000000000.0, 100000000000000000000000.0, 1000000000000000000000000.0, 10000000000000000000000000.0,
            100000000000000000000000000.0, 1000000000000000000000000000.0, 10000000000000000000000000000.0, 100000000000000000000000000000.0, 1000000000000000000000000000000.0,
            10000000000000000000000000000000.0, 100000000000000000000000000000000.0, 1000000000000000000000000000000000.0, 10000000000000000000000000000000000.0, 100000000000000000000000000000000000.0
    };

    while (*s == ' ') s++;
    if (*s == '-') neg = true, s++;
    else if (*s == '+') s++;
    for (sep = -INF; *s; s++) {
        int dig = *s - '0';
        if (dig >= 0 && dig <= 9) {
            tot = 10 * tot + dig, sep++;
            if (tot >= limit) {
                s++;
                if (sep >= 0) goto find_e;
                for (; *s; s++) {
                    if (*s >= '0' && *s <= '9') exponent++;
                    else if (*s == '.' || *s == ',') {
                        s++;
                        goto find_e;
                    } else if (*s == 'e' || *s == 'E') goto at_e;
                    else goto calc;
                }
                goto find_e;
            }
        } else if (*s == '.' || *s == ',') {
            if (sep >= 0) goto calc;
            sep = 0;
        } else if (*s == 'e' || *s == 'E') goto at_e;
        else goto calc;
    }

    find_e:
    for (; *s; s++) {
        if (*s >= '0' && *s <= '9') {} else if (*s == 'e' || *s == 'E') goto at_e;
        else goto calc;
    }
    goto calc;

    at_e:
    s++;
    if (*s) exponent += std::stoi(s);

    calc:
    if (sep > 0) exponent -= sep;
    double res = (neg ? -tot : tot);
    if (exponent > 0) {
        if (exponent <= 35) return res * exps[exponent - 1];
        int num = exponent / 35;
        int extra = (exponent - 35 * num);
        for (int i = 0; i < num; i++) res *= exps[34];
        res *= exps[extra];
        return res;
    } else if (exponent < 0) {
        if (exponent >= -35) return res / exps[-1 - exponent];
        return 0.0;
    }
    return res;
}

inline double to_double(const std::string &s) { return to_double(&(s[0])); }

int unimod_index(const std::string &mod) {
    for (int i = 0; i < UniModIndex.size(); i++) if (mod == UniModIndex[i].first) return UniModIndex[i].second;
    Warning("cannot convert to UniMod: unknown modification");
    return 0;
}

int unimod_index_number(int index) {
    return std::distance(UniModIndices.begin(), std::lower_bound(UniModIndices.begin(), UniModIndices.end(), index));
}

#ifdef WIFFREADER
#ifdef TDFREADER
std::vector<std::string> MSFormatsExt = { ".dia", ".mzml", ".raw", ".tdf", ".wiff" }; // extensions must be ordered alphabetically here
#else
std::vector<std::string> MSFormatsExt = { ".dia", ".mzml", ".raw", ".wiff" };
#endif
HMODULE wiff_dll = NULL;
bool skip_wiff = false, fast_wiff = false;
#else
std::vector<std::string> MSFormatsExt = { ".dia", ".mzml", ".raw" };
#endif

inline double get_mobility(float intensity) { // 1/K0 encoded as part of the intensity 32-bit float
    int enc = *(int*)&intensity;
    enc &= 0xFF;
    float res = enc;
    return (res + 0.5) / 128.0;
}

enum {
    outFile, outRun, outPG, outPID, outPNames, outGenes, outPGQ, outPGN, outGQ, outGN, outGMQ, outGMQU, outModSeq, outStrSeq,
    outPrId, outCharge, outQv, outPQv, outPGQv, outGGQv, outPPt, outPrQ, outPrN, outPrLR, outPrQual,
    outRT, outiRT, outpRT, outpRTf, outpRTt, outpiRT,
    outCols
};

std::vector<std::string> oh = { // output headers
        "File.Name",
        "Run",
        "Protein.Group",
        "Protein.Ids",
        "Protein.Names",
        "Genes",
        "PG.Quantity",
        "PG.Normalised",
        "Genes.Quantity",
        "Genes.Normalised",
        "Genes.MaxLFQ",
        "Genes.MaxLFQ.Unique",
        "Modified.Sequence",
        "Stripped.Sequence",
        "Precursor.Id",
        "Precursor.Charge",
        "Q.Value",
        "Protein.Q.Value",
        "PG.Q.Value",
        "GG.Q.Value",
        "Proteotypic",
        "Precursor.Quantity",
        "Precursor.Normalised",
        "Label.Ratio",
        "Quantity.Quality",
        "RT",
        "RT.Start",
        "RT.Stop",
        "iRT",
        "Predicted.RT",
        "Predicted.iRT",
};

double NormalisationQvalue = 0.01;
double NormalisationPeptidesFraction = 0.4;
bool LocalNormalisation = true;
bool NoRTDepNorm = false;
bool NoSigDepNorm = true;
int LocNormRadius = 250;
double LocNormMax = 2.0;

#define HASH 0
#if (HASH > 0)
unsigned int hash_map[0x10000];
void init_hash() {
	std::mt19937_64 gen(1);
	for (int i = 0; i < 0x10000; i++) hash_map[i] = gen();
}
template<class T> inline unsigned int hashS(T x) { // T mist be 32-bit
	assert(sizeof(T) == 4);
	int y = *((int*)(&x));
	return hash_map[(unsigned int)(y & 0xFFFF)] ^ hash_map[(unsigned int)((y >> 16) & 0xFFFF)];
}
inline unsigned int hashA(int * data, int n) {
	unsigned int res = 0;
	for (int i = 0; i < n; i++) res ^= hash_map[(unsigned int)(data[i] & 0xFFFF)] ^ hash_map[(unsigned int)((data[i] >> 16) & 0xFFFF)];
	return res;
}
unsigned int hashD(float **dataset, int n, int m) {
	unsigned int res = 0;
	for (int i = 0; i < n; i++) res ^= hashA((int*)(dataset[i]), m);
	return res;
}
#endif

bool Predictor = false;
#ifdef PREDICTOR
namespace predictor {
	extern void init_predictor(int id_range, int threads);
	extern bool code_from_precursor(vector<long long> &code, const string &precursor, string &temp_s, vector<string> &temp_sv, const vector<pair<string, int> > &dict);

	class Predictor {
	public:
		Predictor();
		void set_instance(int instance_id);
		map<string, int> get_aa_indices();
		void predict(vector<vector<float> > &spectra, vector<pair<vector<long long>, int> > &codes, int mode, bool verbose = false);
	};
}

predictor::Predictor P;
#endif

class Parameter {
public:
    int min_iter_seek, min_iter_learn;

    Parameter(int _min_iter_seek, int _min_iter_learn) { min_iter_seek = _min_iter_seek, min_iter_learn = _min_iter_learn; }
};

int MaxF = INF, MinF = 0;
const int TopF = 6, auxF = 12;
const int nnS = 2, nnW = (2 * nnS) + 1;
const int qL = 3;
const int QSL[qL] = { 1, 3, 5 };
enum {
    pTimeCorr,
    pLocCorr, pMinCorr, pTotal, pCos, pCosCube, pMs1TimeCorr, pNFCorr, pdRT, pResCorr, pResCorrNorm, pTightCorrOne, pTightCorrTwo, pShadow, pHeavy,
    pMs1TightOne, pMs1TightTwo,
    pMs1Iso, pMs1IsoOne, pMs1IsoTwo,
    pMs1Ratio,
    pBestCorrDelta, pTotCorrSum,
    pMz, pCharge, pLength, pFrNum,
#if AAS
    pMods,
    pAAs,
    pRT = pAAs + 20,
#else
    pRT,
#endif
    pAcc,
    pRef = pAcc + TopF,
    pSig = pRef + auxF - 1,
    pCorr = pSig + TopF,
    pShadowCorr = pCorr + auxF,
    pShape = pShadowCorr + TopF,
#if Q1
    pQLeft = pShape + nnW,
    pQRight = pQLeft + qL,
    pQPos = pQRight + qL,
    pQNFCorr = pQPos + qL,
    pQCorr = pQNFCorr + qL,
    pN = pQCorr + qL
#else
    pN = pShape + nnW
#endif
};

const int TIC_n = 250;
class RunStats {
public:
    RunStats() {}
    int precursors, proteins;
    double tot_quantity, MS1_tot, MS2_tot, MS1_acc, MS1_acc_corrected, MS2_acc, MS2_acc_corrected, RT_acc, pep_length, pep_charge, pep_mc, fwhm_scans, fwhm_RT;
    int cnt_valid;
    float TIC[TIC_n];
};
bool SaveRunStats = true;


std::vector<int> nn_mod, nn_shift, nn_scale;
int nn_mod_n = 4;

std::string location_to_file_name(const std::string& loc) {
    auto result = loc;
    for (int i = 0; i < result.size(); i++) {
        char c = result[i];
        if (c == ':' || c == '/' || c == '\\' || c == '.') result[i] = '_';
    }
    return temp_folder + result;
}

std::string get_extension(const std::string &file) {
    std::string res;
    int pos = file.find_last_of('.');
    if (pos == std::string::npos) return res;
    res = file.substr(pos);
    pos = res.length();
    if (pos > 0) while (res[pos - 1] == '\"') {
            pos--;
            if (!pos) break;
        }
    res = res.substr(0, pos);
    std::transform(res.begin(), res.end(), res.begin(), ::tolower);
    return res;
}

std::string remove_extension(const std::string &file) {
    std::string res;
    int pos = file.find_last_of('.');
    if (pos == std::string::npos) return file;
    res = file.substr(0, pos);
    return res;
}

template <class T> inline double sum(T * x, int n) {
    double ex = 0.0;
    for (int i = 0; i < n; i++) ex += x[i];
    return ex;
}
template <class T> inline double sum(std::vector<T> &x) { return sum(&(x[0]), x.size()); }

template <class T> inline double mean(T * x, int n) {
    assert(n >= 1);
    return sum(x, n) / (double)n;
}
template <class T> inline double mean(std::vector<T> &x) {
    int n = x.size();
    assert(n >= 1);
    return sum(x) / (double)n;
}

template <class T> inline double var(T * x, int n) {
    double ex = mean(x, n), s2 = 0.0;
    for (int i = 0; i < n; i++) s2 += Sqr(x[i] - ex);
    if (n > 1) s2 /= (double)(n - 1);
    return s2;
}

template <class T> inline double var(std::vector<T> &x) { return var(&(x[0]), x.size()); }


template <class Tx, class Ty> inline double cos(Tx * x, Ty * y, int n) {
    double s2x = 0.0, s2y = 0.0, r = 0.0;
    for (int i = 0; i < n; i++) s2x += Sqr(x[i]), s2y += Sqr(y[i]), r += x[i] * y[i];
    if (s2x < E || s2y < E) return 0.0;
    return r / sqrt(s2x * s2y);
}

template <class Tx, class Ty> inline double scalar(Tx * x, Ty * y, int n) {
    double r = 0.0;
    for (int i = 0; i < n; i++) r += x[i] * y[i];
    return r;
}

template <class Tx, class Ty> inline double scalar(Tx * x, Ty * y, bool * mask, int n) {
    double r = 0.0;
    for (int i = 0; i < n; i++) if (mask[i]) r += x[i] * y[i];
    return r;
}

class Group {
public:
    mutable int index, size;
    std::vector<int> e;
    mutable std::vector<int> s;

    Group(std::vector<int> &_e, std::vector<int> &_s) { e = _e, s = _s; }
    inline friend bool operator < (const Group &left, const Group &right) { return left.e < right.e; }
};

std::vector<Group> bipartite_set_cover(std::vector<std::vector<int> > &sets, int N) { // greedy algorithm
    int i, max_size, size, M = sets.size();

    // merge identical sets into groups
    std::set<Group> merged;
    std::vector<Group> groups; // result
    std::vector<int> v(1);
    for (i = 0; i < M; i++) {
        auto &s = sets[i];
        if (!s.size()) continue;
        v[0] = i;
        auto ins = merged.insert(Group(sets[i], v));
        if (!ins.second) ins.first->s.push_back(i);
    }
    if (!merged.size()) return groups;

    // index
    i = max_size = 0;
    int S = merged.size();
    for (auto it = merged.begin(); it != merged.end(); it++) {
        it->index = i++;
        it->size = it->e.size();
        if (it->size > max_size) max_size = it->size;
    }

    // annotate elements
    std::vector<int> removed(N);
    std::vector<std::vector<std::set<Group>::iterator> > elements(N);
    for (auto it = merged.begin(); it != merged.end(); it++) for (auto &e : it->e) elements[e].push_back(it);

    // create size index
    std::vector<int> added(S);
    std::vector<std::list<std::set<Group>::iterator> > size_index(max_size + 1);
    std::vector<std::list<std::set<Group>::iterator>::iterator> index(S);
    for (auto it = merged.begin(); it != merged.end(); it++) {
        auto &s = *it;
        size = s.size;
        size_index[size].push_back(it);
        auto end = size_index[size].end();
        index[it->index] = --end;
    }

    while ((size = size_index.size() - 1) > 0) {
        auto se = size_index.end();
        auto &si = *(--se);
        auto ie = si.end();
        auto &sl = *(--ie); // iterator to one of the groups with the maximum number of elements

        int ind = sl->index;
        added[ind] = 1;
        groups.push_back(*sl);
        si.pop_back();
        auto &e = groups[groups.size() - 1].e;
        for (auto &el : e) if (!removed[el]) {
                removed[el] = 1;
                auto &gr = elements[el];
                for (auto &g : gr) if (!added[g->index]) {
                        size = g->size--;
                        size_index[size].erase(index[g->index]);
                        if (g->size >= 1) size_index[g->size].push_back(g);
                        auto end = size_index[g->size].end();
                        index[g->index] = --end;
                    }
            }
        while (!size_index[size_index.size() - 1].size() && size_index.size()) size_index.resize(size_index.size() - 1);
        if (!size_index.size()) break;
    }

    return groups;
}

class Peak {
public:
    float mz;
    float height;

    Peak() {}
    Peak(float _mz, float _height) {
        mz = _mz;
        height = _height;
    }
    void init(float _mz, float _height) {
        mz = _mz;
        height = _height;
    }
    friend inline bool operator < (const Peak &left, const Peak &right) { return left.mz < right.mz; }

};

struct Trace {
    float mz;
    float height;
    float score;

    Trace() {}
    Trace(float _mz, float _height, float _score) {
        mz = _mz;
        height = _height;
        score = _score;
    }
    void init(float _mz, float _height, float _score) {
        mz = _mz;
        height = _height;
        score = _score;
    }
    friend inline bool operator < (const Trace &left, const Trace &right) { return left.mz < right.mz; }
};

const int fTypeB = 1 << 0;
const int fTypeY = 1 << 1;
const int fExclude = 1 << 6;

class Product {
public:
    float mz = 0.0;
    float height = 0.0;
    char charge = 0, type = 0, index = 0, loss = 0;

    Product() {}
    Product(float _mz, float _height, int _charge) {
        mz = _mz;
        height = _height;
        charge = _charge;
    }
    Product(float _mz, float _height, int _charge, int _type, int _index, int _loss) {
        mz = _mz;
        height = _height;
        charge = _charge;
        type = _type;
        index = _index;
        loss = _loss;
    }
    void init(float _mz, float _height, int _charge) {
        mz = _mz;
        height = _height;
        charge = _charge;
    }
    friend inline bool operator < (const Product &left, const Product &right) { return left.mz < right.mz; }

    inline int ion_code() { return (((((int)type) * 20 + (int)charge)) * (loss_other + 1) + (int)loss) * 100 + (int)index + 1; }
#if (HASH > 0)
    unsigned int hash() { return hashS(mz) ^ hashS(height); }
#endif
};

double AA[256];
int AA_index[256];
const char * AAs = "GAVLIFMPWSCTYHKRQEND";
const char * MutateAAto = "LLLVVLLLLTSSSSLLNDQE";

void init_all() {
    for (int i = 0; i < 256; i++) AA_index[i] = 0, AA[i] = 0.0;

    AA['G'] = 57.021464;
    AA['A'] = 71.037114;
    AA['V'] = 99.068414;
    AA['L'] = 113.084064;
    AA['I'] = 113.084064;

    AA['F'] = 147.068414;
    AA['M'] = 131.040485;
    AA['P'] = 97.052764;
    AA['W'] = 186.079313;
    AA['S'] = 87.032028;

    AA['C'] = 103.009185;
    AA['T'] = 101.047679;
    AA['Y'] = 163.06332;
    AA['H'] = 137.058912;
    AA['K'] = 128.094963;

    AA['R'] = 156.101111;
    AA['Q'] = 128.058578;
    AA['E'] = 129.042593;
    AA['N'] = 114.042927;
    AA['D'] = 115.026943;

    AA['U'] = 150.95363;
    AA['X'] = 196.995499;

    for (int i = 0; i < 20; i++)
        AA_index[AAs[i]] = i;
    AA_index['U'] = AA_index['C'];
    AA_index['X'] = AA_index['M'];
}

inline bool do_cut(char aa, char next) {
    bool cut = false;
    int16_t query = (aa << 8) | '*';
    for (auto &s : Cut) if (s == query) { cut = true; break; }
    if (!cut) {
        query = ('*' << 8) | next;
        for (auto &s : Cut) if (s == query) { cut = true; break; }
        if (!cut) {
            query = (aa << 8) | next;
            for (auto &s : Cut) if (s == query) { cut = true; break; }
        }
    }
    if (cut && NoCut.size()) {
        query = ('*' << 8) | next;
        for (auto &s : NoCut) if (s == query) { cut = false; break; }
        if (cut) {
            query = ('*' << 8) | next;
            for (auto &s : NoCut) if (s == query) { cut = false; break; }
            if (cut) {
                query = (aa << 8) | next;
                for (auto &s : NoCut) if (s == query) { cut = false; break; }
            }
        }
    }
    return cut;
}

std::vector<double> yDelta;
int yDeltaS = 2, yDeltaW = (yDeltaS * 2 + 1);
inline int y_delta_composition_index(char aa) { return AA_index[aa]; }
inline int y_delta_charge_index() { return y_delta_composition_index(AAs[19]) + 1; }
inline int y_delta_index(char aa, int shift) {
    assert(Abs(shift) <= yDeltaS);
    return y_delta_charge_index() + AA_index[aa] * yDeltaW + shift + yDeltaS + 1;
}
int yCTermD = 9, yNTermD = 8;
inline int y_cterm_index(char aa, int shift) {
    int row = (aa == 'K' ? 0 : (aa == 'R' ? 1 : 2));
    return row * yCTermD + Min(shift, yCTermD - 1) + y_delta_index(AAs[19], yDeltaS) + 1;
}
inline int y_nterm_index(int shift) { return Min(shift, yNTermD - 1) + y_cterm_index(AAs[19], yCTermD - 1) + 1; }
inline int y_delta_size() { return y_nterm_index(yNTermD - 1) + 1; }

inline double y_ratio(int i, int charge, const std::string &aas) {
    int n = aas.length(), j;
    double v = yDelta[y_delta_charge_index()] * (double)charge;
    for (j = 0; j < n; j++) v += yDelta[y_delta_composition_index(aas[j])];
    for (j = Max(i - yDeltaS, 0); j <= Min(i + yDeltaS, n - 1); j++)
        v += yDelta[y_delta_index(aas[j], j - i)];
    v += yDelta[y_nterm_index(i - 2)];
    v += yDelta[y_cterm_index(aas[n - 1], n - 1 - i)];

    return v;
}

void y_scores(std::vector<double> &s, int charge, const std::string &aas) {
    int n = aas.length(), i;
    s.resize(n);
    s[0] = s[1] = 0.0;
    for (i = 2; i < n; i++)
        s[i] = s[i - 1] + y_ratio(i, charge, aas);
}

void to_exp(std::vector<double> &s) { for (int i = 0; i < s.size(); i++) s[i] = exp(s[i]); }

std::vector<double> yLoss;
inline int y_loss_aa(char aa, bool NH3) { return AA_index[aa] + 20 * (int)NH3; }
inline int y_loss(bool NH3) { return y_loss_aa(AAs[19], true) + 1 + (int)NH3; }
inline int y_loss_size() { return y_loss(true) + 1; }

void y_loss_scores(std::vector<double> &s, const std::string &aas, bool NH3) {
    int n = aas.length(), i;
    double sum = yLoss[y_loss(NH3)];
    s.resize(n); s[0] = 0.0;
    for (i = n - 1; i >= 1; i--) {
        sum += yLoss[y_loss_aa(aas[i], NH3)];
        s[i] = sum;
    }
}

std::vector<double> InSilicoRT;
int RTTermD = 1, RTTermDScaled = 5;
inline int aa_rt_nterm(char aa, int shift) { return 1 + Min(shift, RTTermD - 1) * 20 + AA_index[aa]; }
inline int aa_pointsterm(char aa, int shift) { return aa_rt_nterm(AAs[19], RTTermD - 1) + 1 + Min(shift, RTTermD - 1) * 20 + AA_index[aa]; }
inline int aa_rt_nterm_scaled(char aa, int shift) { return aa_pointsterm(AAs[19], RTTermD - 1) + 1 + Min(shift, RTTermDScaled - 1) * 20 + AA_index[aa]; }
inline int aa_pointsterm_scaled(char aa, int shift) { return aa_rt_nterm_scaled(AAs[19], RTTermDScaled - 1) + 1 + Min(shift, RTTermDScaled - 1) * 20 + AA_index[aa]; }
inline int mod_rt_index(int mod) { return aa_pointsterm_scaled(AAs[19], RTTermDScaled - 1) + 1 + unimod_index_number(mod); }
inline int mod_rt_index(const std::string &mod) { return mod_rt_index(unimod_index(mod)); }
inline int in_silico_rt_size() { return UniModIndices.size() ? (mod_rt_index(UniModIndices[UniModIndices.size() - 1]) + 1) : (aa_pointsterm_scaled(AAs[19], RTTermDScaled - 1) + 1); }

void init_prediction() {
    yDelta.resize(y_delta_size(), 0.0);
    yLoss.resize(y_loss_size(), 0.0);
    InSilicoRT.resize(in_silico_rt_size(), 0.0);
}

inline int closing_bracket(const std::string &name, char symbol, int pos) {
    int end, par;
    char close = (symbol == '(' ? ')' : ']');
    for (end = pos + 1, par = 1; end < name.size(); end++) {
        char s = name[end];
        if (s == close) {
            par--;
            if (!par) break;
        } else if (s == symbol) par++;
    }
    return end;
}

inline int peptide_length(const std::string &name) {
    int i, n = 0;
    for (i = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i = closing_bracket(name, symbol, i);
            continue;
        }
        n++;
    }
    return n;
}

inline int missed_KR_cleavages(const std::string &name) {
    int i, n = 0, last = 0, par, end;
    for (i = 0; i < name.size() - 1; i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i = closing_bracket(name, symbol, i);
            continue;
        }
        if (symbol == 'K' || symbol == 'R') n++, last = 1;
        else last = 0;
    }
    return n - last;
}

inline std::string get_aas(const std::string &name) {
    int i, end, par;
    std::string result;

    for (i = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i = closing_bracket(name, symbol, i);
            continue;
        }
        result.push_back(symbol);
    }
    return result;
}

inline std::vector<bool> get_mod_state(const std::string &name, int length) {
    int i, j;
    std::vector<bool> result(length, false);

    for (i = j = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            result[j] = true;
            i = closing_bracket(name, symbol, i);
            continue;
        }
        j++;
    }
    return result;
}

inline std::vector<double> get_sequence(const std::string &name, int * no_cal = NULL) {
    int i, j;
    double add = 0.0;
    std::vector<double> result;
    if (no_cal != NULL) *no_cal = 0;

    for (i = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i++;

            int end = closing_bracket(name, symbol, i);
            if (end >= name.size()) Throw(QStringLiteral("incorrect peptide name format: " ) );

            std::string mod = name.substr(i, end - i);
            for (j = 0; j < Modifications.size(); j++) if (Modifications[j].name == mod) {
                    if (!result.size()) add += Modifications[j].mass;
                    else result[result.size() - 1] += Modifications[j].mass;
                    if (no_cal != NULL && Modifications[j].mass < 4.5 && Modifications[j].mass > 0.0 && Abs(Modifications[j].mass - (double)floor(Modifications[j].mass + 0.5)) < 0.1) *no_cal = 1;
                    break;
                }
            if (j == Modifications.size())
            qDebug() << QStringLiteral("unknown modification: ");
            i = end;
            continue;
        }

        result.push_back(AA[symbol] + add);
        add = 0.0;
    }
    return result;
}

inline std::string to_canonical(const std::string &name) {
    int i, j;
    std::string result;

    for (i = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i++;

            int end = closing_bracket(name, symbol, i);
            if (end == std::string::npos) qDebug() << QStringLiteral("incorrect peptide name format: ");

            std::string mod = name.substr(i, end - i);
            for (j = 0; j < UniMod.size(); j++) if (UniMod[j].first == mod) {
                    result += std::string("(") + UniMod[j].second + std::string(")");
                    break;
                }
            if (j == Modifications.size()) result += std::string("(") + mod + std::string(")"); // no warning here
            i = end;
            continue;
        } else result.push_back(symbol);
    }
    return result;
}
inline std::string to_canonical(const std::string &name, int charge) { return to_canonical(name) + std::to_string(charge); }

inline void to_eg(std::string &eg, const std::string &name) {
    int i, j;

    eg.clear();
    for (i = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i++;

            int end = closing_bracket(name, symbol, i);
            if (end == std::string::npos) qDebug() << QStringLiteral("incorrect peptide name format: ");

            std::string mod = name.substr(i, end - i);
            for (j = 0; j < Modifications.size(); j++) if (Modifications[j].name == mod) {
                    if (!Modifications[j].label) eg += std::string("(") + UniMod[j].second + std::string(")");
                    break;
                }
            if (j == Modifications.size()) eg += std::string("(") + mod + std::string(")"); // no warning here
            i = end;
            continue;
        } else eg.push_back(symbol);
    }
}

inline std::string to_eg(const std::string &name) {
    std::string eg;
    to_eg(eg, name);
    return eg;
}

inline std::string to_charged_eg(const std::string &name, int charge) { return to_eg(name) + std::to_string(charge); }

inline std::string to_prosit(const std::string &name) {
    int i, j;
    std::string result;

    for (i = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i++;

            int end = closing_bracket(name, symbol, i);
            if (end == std::string::npos) qDebug() << QStringLiteral("incorrect peptide name format: ");

            if (i + 9 < name.size()) if (!memcmp(&(name[i]), "UniMod:35", 9)) result += "(ox)"; // methionine oxidation
            i = end;
            continue;
        } else result.push_back(symbol);
    }
    return result;
}

inline std::string pep_name(const std::string &pr_name) {
    int i;
    std::string result = pr_name;
    for (i = result.size() - 1; i >= 0; i--) if (result[i] < '0' || result[i] > '9') break;
    result.resize(i + 1);
    return result;
}

void count_rt_aas(std::vector<double> &count, const std::string &name) {
    int i, pos, l = peptide_length(name);
    double scale = 1.0 / (double)l;

    for (i = pos = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i++;

            int end = closing_bracket(name, symbol, i);
            if (end == std::string::npos) qDebug() << QStringLiteral("incorrect peptide name format: ");
            count[mod_rt_index(name.substr(i, end - i))]++;
            i = end;
            continue;
        }
        count[aa_rt_nterm(symbol, pos)] += 1.0;
        count[aa_pointsterm(symbol, l - 1 - pos)] += 1.0;
        count[aa_rt_nterm_scaled(symbol, pos)] += scale;
        count[aa_pointsterm_scaled(symbol, l - 1 - pos)] += scale;
        pos++;
    }
}

double predict_irt(const std::string &name) {
    int i, pos, l = peptide_length(name);
    double iRT = InSilicoRT[0], scale = 1.0 / (double)l;

    for (i = pos = 0; i < name.size(); i++) {
        char symbol = name[i];
        if (symbol < 'A' || symbol > 'Z') {
            if (symbol != '(' && symbol != '[') continue;
            i++;

            int end = closing_bracket(name, symbol, i);
            if (end == std::string::npos) qDebug() << QStringLiteral("incorrect peptide name format: ");
            iRT += InSilicoRT[mod_rt_index(name.substr(i, end - i))];
            i = end;
            continue;
        }
        iRT += InSilicoRT[aa_rt_nterm(symbol, pos)];
        iRT += InSilicoRT[aa_pointsterm(symbol, l - 1 - pos)];
        iRT += InSilicoRT[aa_rt_nterm_scaled(symbol, pos)] * scale;
        iRT += InSilicoRT[aa_pointsterm_scaled(symbol, l - 1 - pos)] * scale;
        pos++;
    }
    return iRT;
}

template <class F, class T> void write_vector(F &out, std::vector<T> &v) {
    int size = v.size();
    out.write((char*)&size, sizeof(int));
    if (size) out.write((char*)&(v[0]), size * sizeof(T));
}

template <class F, class T> void read_vector(F &in, std::vector<T> &v) {
    int size = 0; in.read((char*)&size, sizeof(int));
    if (size) {
        v.resize(size);
        in.read((char*)&(v[0]), size * sizeof(T));
    }
}

template<class F> void write_string(F &out, std::string &s) {
    int size = s.size();
    out.write((char*)&size, sizeof(int));
    out.write((char*)&(s[0]), size);
}

template<class F> void read_string(F &in, std::string &s) {
    int size = 0; in.read((char*)&size, sizeof(int));
    if (size) {
        s.resize(size);
        in.read((char*)&(s[0]), size);
    }
}

template <class F, class T> void write_array(F &out, std::vector<T> &a) {
    int size = a.size();
    out.write((char*)&size, sizeof(int));
    for (int i = 0; i < size; i++) a[i].write(out);
}

template <class F, class T> void read_array(F &in, std::vector<T> &a) {
    int size = 0; in.read((char*)&size, sizeof(int));
    if (size) {
        a.resize(size);
        for (int i = 0; i < size; i++) a[i].read(in);
    }
}

template<class F> void write_strings(F &out, std::vector<std::string> &strs) {
    int size = strs.size();
    out.write((char*)&size, sizeof(int));
    for (int i = 0; i < size; i++) write_string(out, strs[i]);
}

template<class F> void read_strings(F &in, std::vector<std::string> &strs) {
    int size = 0; in.read((char*)&size, sizeof(int));
    if (size) {
        strs.resize(size);
        for (int i = 0; i < size; i++) read_string(in, strs[i]);
    }
}

char char_from_type[3] = { '?', 'b', 'y' };
std::string name_from_loss[loss_other + 1] = { std::string("noloss"), std::string("H2O"), std::string("NH3"), std::string("CO"), std::string("unknown"), std::string("unknown") };

class Ion {
public:
    char type, index, loss, charge; // index - number of AAs to the left from the cleavage site
    float mz;

    Ion() {};
    Ion(int _type, int _index, int _loss, int _charge, double _mz) {
        type = _type;
        index = _index;
        loss = _loss;
        charge = _charge;
        mz = _mz;
    }

    template<class T> void init(T &other) {
        type = other.type & 3;
        index = other.index;
        loss = other.loss;
        charge = other.charge;
        mz = other.mz;
    }

    bool operator == (const Ion &other) { return index == other.index && type == other.type && charge == other.charge && loss == other.loss; }
};

std::vector<Ion> generate_fragments(const std::vector<double> &sequence, int charge, int loss, int *cnt, int min_aas = 0) {
    int i;
    double c = (double)charge, curr, s = sum(&(sequence[0]), sequence.size());
    std::vector<Ion> result;

    for (i = 0, curr = 0.0; i < sequence.size() - 1; i++) {
        curr += sequence[i];
        double b = curr;
        double y = s - curr + proton + OH;

        if (i + 1 >= min_aas) result.push_back(Ion(fTypeB, i + 1, loss, charge, (b + c * proton - Loss[loss]) / c)), (*cnt)++;
        if (sequence.size() - i - 1 >= min_aas) result.push_back(Ion(fTypeY, i + 1, loss, charge, (y + c * proton - Loss[loss]) / c)), (*cnt)++;
    }
    return result;
}

std::vector<Ion> generate_all_fragments(const std::vector<double> &sequence, int charge, int max_loss, int *cnt) {
    int i;
    double c = (double)charge, curr, s = sum(&(sequence[0]), sequence.size());
    std::vector<Ion> result;

    for (int loss = loss_none; loss <= max_loss; loss++) {
        for (i = 0, curr = 0.0; i < sequence.size() - 1; i++) {
            curr += sequence[i];
            double b = curr;
            double y = s - curr + proton + OH;

            result.push_back(Ion(fTypeB, i + 1, loss, charge, (b + c * proton - Loss[loss]) / c)), (*cnt)++;
            result.push_back(Ion(fTypeY, i + 1, loss, charge, (y + c * proton - Loss[loss]) / c)), (*cnt)++;
        }
    }
    return result;
}

double max_gen_acc = 0.0;
Lock gen_acc_lock;
std::vector<Ion> recognise_fragments(const std::vector<double> &sequence, const std::vector<Product> &fragments, float &gen_acc, int pr_charge, bool full_spectrum = false, int max_charge = 19, int loss_cap = loss_N) {
    int i, j, cnt, tot, charge, loss, index;
    std::vector<Ion> result(fragments.size());
    std::vector<Ion> v;
    double delta, min, s = sum(&(sequence[0]), sequence.size());
    gen_acc = GeneratorAccuracy;

    anew:
    cnt = tot = index = 0, charge = 1;
    loss = loss_none;
    for (i = 0; i < result.size(); i++) result[i].charge = 0;

    start:
    v = generate_fragments(sequence, charge, loss, &cnt);
    if (charge == 1) v.push_back(Ion(fTypeY, sequence.size(), loss, pr_charge, (s + OH + (1.0 + double(pr_charge)) * proton - Loss[loss]) / double(pr_charge))); // non-fragmented

    for (i = 0; i < fragments.size(); i++) {
        if (result[i].charge) continue;
        double margin = fragments[i].mz * gen_acc;
        for (j = 0, min = margin; j < v.size(); j++) {
            delta = Abs(v[j].mz - fragments[i].mz);
            if (delta < min) min = delta, index = j;
        }
        if (min < margin) {
            tot++;
            result[i].init(v[index]);
            if (tot >= fragments.size()) goto stop;
        }
    }

    loss++;
    if (loss == loss_cap) loss = loss_none, charge++;
    if (charge <= max_charge && tot < fragments.size()) goto start;

    stop:
    if (tot < fragments.size()) {
        if (full_spectrum) return result;
        gen_acc *= 2.0;
        if (Verbose >= 1 && gen_acc > max_gen_acc)
            qDebug() << "WARNING: not all fragments recognised; generator accuracy threshold increased to " << gen_acc << "\n";
        if (gen_acc > max_gen_acc) {
            while (!gen_acc_lock.set()) {}
            max_gen_acc = gen_acc;
            gen_acc_lock.free();
        }
        v.clear();
        goto anew;
    }

    return result;
}

std::vector<Ion> recognise_fragments(const std::string &name, const std::vector<Product> &fragments, float &gen_acc, int pr_charge, bool full_spectrum = false, int max_charge = 19, int loss_cap = loss_N) {
    int i;
    if (fragments.size()) if (fragments[0].type & 3) {
            std::vector<Ion> result(fragments.size());
            for (i = 0; i < result.size(); i++)
                if (fragments[i].type & 3) result[i].init(fragments[i]); else break;
            if (i == result.size()) return result;
        }
    auto sequence = get_sequence(name);
    return recognise_fragments(sequence, fragments, gen_acc, pr_charge, full_spectrum, max_charge, loss_cap);
}


std::vector<Ion> generate_fragments(const std::vector<double> &sequence, const std::vector<Ion> &pattern, int pr_charge) {
    int i, j, charge = 1, loss = loss_none, cnt = 0, tot = 0;
    double s = sum(&(sequence[0]), sequence.size());
    std::vector<Ion> result(pattern.size());
    std::vector<Ion> v;

    start:
    v = generate_fragments(sequence, charge, loss, &cnt); i = 0;
    if (charge == 1) v.push_back(Ion(fTypeY, sequence.size(), loss, pr_charge, (s + OH + (1.0 + double(pr_charge)) * proton - Loss[loss]) / double(pr_charge))); // non-fragmented

    for (auto pos = pattern.begin(); pos != pattern.end(); pos++, i++) {
        for (j = 0; j < v.size(); j++) if (v[j] == *pos) {
                result[i].init(v[j]); tot++;
                break;
            }
    }
    loss++;
    if (loss == loss_N) loss = loss_none, charge++;
    if (charge < 20 && tot < pattern.size()) goto start;

    return result;
}

bool ProfileSpectrumWarning = false;

struct FI {
    float value;
    int index;

    friend bool inline operator < (const FI& left, const FI& right) { return left.value < right.value; }
    friend bool inline operator > (const FI& left, const FI& right) { return left.value > right.value; }
};

struct Mass {
    float mz;
    int pr, index, charge;

    Mass() {}
    Mass(float _mz, int _pr, int _index, int _charge) { mz = _mz, pr = _pr, index = _index, charge = _charge; }
    friend inline bool operator < (const Mass& left, const Mass& right) { return left.mz < right.mz; }
};

struct Hit {
    float score;
    int fr;
};

const int SCSPeaks = 50;

const int SCSFr = 12;
struct Match {
    Hit hit[SCSFr];

    inline double score() { int i; double s; for (i = 0, s = 0.0; i < SCSFr; i++) s += hit[i].score; return s; }
};

struct Feature {
    double mz = 0.0, RT = 0.0, score = 0.0, height = 0.0, quantity = 0.0, qvalue = 0.0, search[2], mod[2];
    int apex = 0, ms1_apex = 0, charge = 0, shadow = 0, index[2];

    std::vector<Trace> peaks;

    template <class F> void write(F &out, int scan_number) {
        out << "BEGIN IONS\nTITLE=ssc." << scan_number << '.' << scan_number << ".0\nPEPMASS=" << (float)mz << ' ' << quantity << "\nRTINSECONDS=" << (float)(RT * 60.0) << '\n';
        for (int i = 0; i < peaks.size(); i++) out << peaks[i].mz << ' ' << peaks[i].height << '\n';
        out << "END IONS\n";
    }

    friend inline bool operator < (const Feature &left, const Feature &right) {
        if (left.RT < right.RT) return true;
        if (left.RT <= right.RT && left.mz < right.mz) return true;
        return false;
    }

    template <bool open> void find(bool _decoy, std::vector<Mass> &fragments, std::vector<Mass> &fragments_shift,
                                   std::vector<std::pair<float, int> > &precursors, double acc_ms1, double acc_ms2, std::vector<Match> &temp) {

        temp.clear();
        int decoy = _decoy, frs = fragments.size();
        if (!frs || shadow || charge < 1 || charge > 7) {
            search[decoy] = -INF, index[decoy] = -1, mod[decoy] = 0.0;
            return;
        }

        int i, size = 0, ind;
        float min = mz - mz * acc_ms1, max = mz + mz * acc_ms1, s, sc;
        Mass query(0.0, 0, 0, 0);

        for (int iter = 0; iter <= (int)open; iter++) for (auto &p : peaks) {
                float qmz = ((iter == 0 || !open) ? p.mz : (mz - proton) * ((double)charge) - (p.mz - proton));
                auto &frm = ((iter == 0 || !open) ? fragments : fragments_shift);

                float delta = qmz * acc_ms2;
                query.mz = qmz - delta;
                float high = qmz + delta;
                if (high <= frm[0].mz || query.mz >= frm[frs - 1].mz) continue;

                auto pos = std::lower_bound(frm.begin(), frm.end(), query);
                if (pos == frm.end()) continue;
                while (pos->mz < high) {
                    if (open && iter == 1 && pos->charge != 1) {
                        pos++;
                        if (pos == frm.end()) break;
                    }
                    auto &pr = precursors[pos->pr];
                    if ((open || (pr.first > min && pr.first < max)) && pr.second == charge) {
                        if (size <= pos->pr) size = ((pos->pr * 5) / 4) + 1, temp.resize(size);
                        auto &match = temp[pos->pr].hit;
                        if (p.score <= match[SCSFr - 1].score) goto finish;
                        int fr, j; if (open) fr = pos->index;

                        for (i = 0; i < SCSFr - 1; i++) {
                            auto &h = match[i];
                            if (h.score < E || (open && h.fr == fr && p.score > h.score)) {
                                for (j = i; j > 0; j--) {
                                    auto &next = match[j - 1];
                                    if (next.score >= p.score) break;
                                    match[j] = next;
                                }
                                match[j].score = p.score;
                                if (open) match[j].fr = fr;
                                break;
                            } else if (open && h.fr == fr && p.score <= h.score) break;
                        }
                        if (i == SCSFr - 1) {
                            for (j = i; j > 0; j--) {
                                auto &next = match[j - 1];
                                if (next.score >= p.score) break;
                                match[j] = next;
                            }
                            match[j].score = p.score;
                            if (open) match[j].fr = fr;
                        }
                    }
                    finish:
                    pos++;
                    if (pos == frm.end()) break;
                }
            }

        s = E, ind = -1;
        for (i = 0; i < size; i++) if ((sc = temp[i].score()) > s)
                s = sc, ind = i;
        search[decoy] = s;
        index[decoy] = ind;
        mod[decoy] = (mz - precursors[ind].first) * (double)precursors[ind].second;
    }
};

struct Scan {
    double RT = 0.0, window_low = 0.0, window_high = 0.0;
    int n = 0, type = 0;
    Peak * peaks = NULL;

    Scan() { }

    inline int size() { return n; }

    inline friend bool operator < (const Scan& left, const Scan& right) { return left.RT < right.RT; }

    inline bool has(float mz) { return (window_low <= mz && window_high > mz); }

    inline bool has(float mz, float margin) { return (window_low <= mz + margin && window_high > mz - margin); }

    template <bool get_mz> inline double level(int &low, int &high, float mz, float accuracy, float * peak_mz = NULL) {
        int i;
        float v, s, margin, min, max;
        margin = mz * accuracy;
        min = mz - margin, max = mz + margin;
        if (!size()) {
            if (get_mz) *peak_mz = 0.0;
            return 0.0;
        }

        while (high > low) {
            int middle = (high + low) >> 1;
            v = peaks[middle].mz;
            if (v < max) {
                if (v > min) {
                    s = peaks[middle].height;
                    if (get_mz) *peak_mz = peaks[middle].mz;
                    for (i = middle + 1; i < high && peaks[i].mz < max; i++) if (peaks[i].height > s) {
                            s = peaks[i].height;
                            if (get_mz) *peak_mz = peaks[i].mz;
                        }
                    if (i < high) high = i;
                    for (i = middle - 1; i >= low && peaks[i].mz > min; i--) if (peaks[i].height > s) {
                            s = peaks[i].height;
                            if (get_mz) *peak_mz = peaks[i].mz;
                        }
                    if (i >= low) low = i + 1;
                    return s;
                }
                low = middle + 1;
            } else high = middle;
        }
        return 0.0;
    }

    template <bool get_mz> inline double level(float mz, float accuracy, float * peak_mz = NULL) {
        int low = 0, high = size();
        return level<get_mz>(low, high, mz, accuracy, peak_mz);
    }

    inline double level(float mz, float accuracy, float ref_height) { // finds peak with height closest to the ref_height
        int i, low = 0, high = size();
        float v, s, delta, u, margin, min, max;
        margin = mz * accuracy;
        min = mz - margin, max = mz + margin;
        if (!size()) return 0.0;

        while (high > low) {
            int middle = (high + low) >> 1;
            v = peaks[middle].mz;
            if (v < max) {
                if (v > min) {
                    s = peaks[middle].height;
                    delta = Abs(s - ref_height);
                    for (i = middle + 1; i < high && peaks[i].mz < max; i++) if ((u = Abs(peaks[i].height - ref_height)) < delta) {
                            s = peaks[i].height;
                            delta = u;
                        }
                    if (i < high) high = i;
                    for (i = middle - 1; i >= low && peaks[i].mz > min; i--) if ((u = Abs(peaks[i].height - ref_height)) < delta) {
                            s = peaks[i].height;
                            delta = u;
                        }
                    if (i >= low) low = i + 1;
                    return s;
                }
                low = middle + 1;
            } else high = middle;
        }
        return 0.0;
    }

#if (HASH > 0)
    unsigned int hash() {
		unsigned int res = hashS((float)RT) ^ hashS((float)window_low) ^ hashS((float)window_high);
		for (int i = 0; i < n; i++) res ^= peaks[i].hash();
		return res;
	}
#endif

    void bin_peaks(FI * binned, int bins_per_Da) {
        int i;
        FI * fi = NULL;
        float bpd = (float)bins_per_Da, bpdi = 1.0 / bpd, max = -INF;

        if (!n) return;
        for (i = 0; i < n; i++) {
            auto &p = peaks[i];
            if (p.mz >= max) {
                int bin = floor(p.mz * bpd);
                max = bpdi * float(bin + 1);
                fi = binned + bin;
            }
            if (p.height > fi->value) fi->value = p.height, fi->index = i;
        }
    }

    void bin_peaks(std::vector<FI>& binned, int bins_per_Da) {
        if (!n) return;
        binned.resize(floor(peaks[n - 1].mz * (float)bins_per_Da) + 1);
        bin_peaks(&(binned[0]), bins_per_Da);
    }

    void top_peaks(Peak * top, int N, std::vector<Peak>& temp, float DR = 10000.0) {
        int k;
        float margin;

        if (!n) {
            k = 0;
            goto finish;
        }

        temp.resize(n);
        memcpy(&(temp[0]), peaks, n * sizeof(Peak));
        std::sort(temp.begin(), temp.end(), [&](const Peak &left, const Peak &right) { return left.height > right.height; });

        margin = temp[0].height / DR;
        for (k = 0; k < Min(n, N); k++) if (temp[k].height < margin) break;

        temp.resize(k);
        std::sort(temp.begin(), temp.end());
        memcpy(top, &(temp[0]), k * sizeof(Peak));

        finish:
        if (N > k) memset(top + k, 0, (N - k) * sizeof(Peak));
    }

    void top_peaks(std::vector<Peak> &top, int N, std::vector<Peak>& temp, float DR = 10000.0) {
        top.resize(N);
        top_peaks(&(top[0]), N, temp, DR);
    }

    void top_peaks(Peak * top, int N, float low, float high, std::vector<Peak>& temp, float DR = 10000.0) {
        int i, k;
        float margin;

        if (!n) {
            k = 0;
            goto finish;
        }

        temp.resize(n);
        memcpy(&(temp[0]), peaks, n * sizeof(Peak));
        std::sort(temp.begin(), temp.end(), [&](const Peak &left, const Peak &right) { return left.height > right.height; });

        margin = temp[0].height / DR;
        for (i = k = 0; i < n && k < N; i++) if (temp[i].mz > low && temp[i].mz < high) {
                if (temp[k].height < margin) break;
                top[k++] = temp[i];
            }

        temp.resize(k);
        memcpy(&(temp[0]), top, k * sizeof(Peak));
        std::sort(temp.begin(), temp.end());
        memcpy(top, &(temp[0]), k * sizeof(Peak));

        finish:
        if (N > k) memset(top + k, 0, (N - k) * sizeof(Peak));
    }

    void top_peaks(std::vector<Peak> &top, int N, float low, float high, std::vector<Peak>& temp, float DR = 10000.0) {
        top.resize(N);
        top_peaks(&(top[0]), N, low, high, temp, DR);
    }
};

inline int mass_bin(float mz, double acc) { return floor(log(mz) / acc); }
inline float bin_to_mass(int bin, double acc) { return exp(acc * (double)bin); }
inline float binned_float(float mz, double acc) { return bin_to_mass(mass_bin(mz, acc), acc); }

class Isoform {
public:
    std::string id;
    mutable std::string name, gene, description;
    mutable std::set<int> precursors; // precursor indices in the library
    int name_index = 0, gene_index = 0;
    bool swissprot = true;

    Isoform() {}
    Isoform(const std::string &_id) { id = _id; }
    Isoform(const std::string &_id, const std::string &_name, const std::string &_gene, const std::string &_description, bool _swissprot) {
        id = _id;
        name = _name;
        gene = _gene;
        description = _description;
        swissprot = _swissprot;
    }

    friend inline bool operator < (const Isoform &left, const Isoform &right) { return left.id < right.id; }

    template <class F> void write(F &out) {
        int sp = swissprot, size = precursors.size();
        out.write((char*)&sp, sizeof(int));
        out.write((char*)&size, sizeof(int));
        write_string(out, id);
        write_string(out, name);
        write_string(out, gene);
        out.write((char*)&name_index, sizeof(int));
        out.write((char*)&gene_index, sizeof(int));
        for (auto it = precursors.begin(); it != precursors.end(); it++) out.write((char*)&(*it), sizeof(int));
    }

    template <class F> void read(F &in) {
        int sp = 0, size = 0;
        in.read((char*)&sp, sizeof(int));
        in.read((char*)&size, sizeof(int));
        swissprot = sp;

        read_string(in, id);
        read_string(in, name);
        read_string(in, gene);
        in.read((char*)&name_index, sizeof(int));
        in.read((char*)&gene_index, sizeof(int));
        precursors.clear();
        for (int i = 0; i < size; i++) {
            int pr = -1;
            in.read((char*)&pr, sizeof(int));
            if (pr >= 0) precursors.insert(pr);
        }
    }
};

class PG {
public:
    std::string ids;
    mutable std::string names, genes;
    mutable std::vector<int> precursors; // precursor indices in the library
    mutable std::set<int> proteins;
    std::vector<int> name_indices, gene_indices;

    PG() {}

    PG(std::string _ids) {
        ids = _ids;
    }

    friend inline bool operator < (const PG &left, const PG &right) { return left.ids < right.ids; }

    void annotate(std::vector<Isoform> &_proteins, std::vector<std::string> &_names, std::vector<std::string> &_genes) {
        std::set<int> name, gene;
        std::string word;
        if (_names.size()) {
            for (auto &p : proteins) if (_names[_proteins[p].name_index].size()) name.insert(_proteins[p].name_index);
            if (names.size() && !proteins.size()) {
                std::stringstream list(names);
                while (std::getline(list, word, ';')) {
                    word = fast_trim(word);
                    auto pos = std::lower_bound(_names.begin(), _names.end(), word);
                    if (pos != _names.end()) if (*pos == word) name.insert(std::distance(_names.begin(), pos));
                }
            }
            names.clear(); name_indices.clear(); name_indices.insert(name_indices.begin(), name.begin(), name.end());
            if (name_indices.size()) for (int i = 0; i < name_indices.size(); i++) names += (names.size() ? std::string(";") : std::string("")) + _names[name_indices[i]];
        }
        if (_genes.size()) {
            for (auto &p : proteins) if (_genes[_proteins[p].gene_index].size()) gene.insert(_proteins[p].gene_index);
            if (genes.size() && !proteins.size()) {
                std::stringstream list(genes);
                while (std::getline(list, word, ';')) {
                    word = fast_trim(word);
                    auto pos = std::lower_bound(_genes.begin(), _genes.end(), word);
                    if (pos != _genes.end()) if (*pos == word) gene.insert(std::distance(_genes.begin(), pos));
                }
            }
            genes.clear(); gene_indices.clear(); gene_indices.insert(gene_indices.begin(), gene.begin(), gene.end());
            if (gene_indices.size()) for (int i = 0; i < gene_indices.size(); i++) genes += (genes.size() ? std::string(";") : std::string("")) + _genes[gene_indices[i]];
        }
    }

    template <class F> void write(F &out) {
        int size_p = proteins.size();
        out.write((char*)&size_p, sizeof(int));
        write_string(out, ids);
        write_string(out, names);
        write_string(out, genes);
        write_vector(out, name_indices);
        write_vector(out, gene_indices);
        write_vector(out, precursors);
        for (auto it = proteins.begin(); it != proteins.end(); it++) out.write((char*)&(*it), sizeof(int));
    }

    template <class F> void read(F &in) {
        int sp = 0, size_p = 0;
        in.read((char*)&size_p, sizeof(int));

        read_string(in, ids);
        read_string(in, names);
        read_string(in, genes);
        read_vector(in, name_indices);
        read_vector(in, gene_indices);
        read_vector(in, precursors);
        for (int i = 0; i < size_p; i++) {
            int p = -1;
            in.read((char*)&p, sizeof(int));
            if (p >= 0) proteins.insert(p);
        }
    }
};

enum {
    qTotal, qFiltered,
    qN
};

class Fragment {
public:
    mutable int index; // index of the fragment in the library (i.e. its number among the fragments of the precursor)
    mutable float quantity[qN];
    mutable float corr;
};

class PrecursorEntry {
public:
    int index; // index of the precursor in the spectral library
    int run_index;
    mutable bool decoy_found;
    mutable int apex, peak, best_fragment, peak_width;
    mutable float RT, RT_start, RT_stop, iRT, predicted_RT, predicted_iRT, one_over_K0, best_fr_mz, profile_qvalue, qvalue, lib_qvalue, dratio, decoy_qvalue, protein_qvalue, pg_qvalue, gg_qvalue, quantity, ratio, level, norm, quality;
    mutable float pg_quantity, pg_norm, gene_quantity, gene_norm, max_lfq, max_lfq_unique;
    mutable float evidence, decoy_evidence, ms1_corr, ms1_area, cscore, decoy_cscore;
#if REPORT_SCORES
    float scores[pN], decoy_scores[pN];
#endif
    friend inline bool operator < (const PrecursorEntry &left, const PrecursorEntry &right) { return left.index < right.index; }
};

class QuantEntry {
public:
    int index = -1, window, fr_n = 0; // index must be initialised with -1
    mutable PrecursorEntry pr;
    mutable Fragment fr[TopF];
#if ELUTION_PROFILE
    mutable std::pair<float, float> ms1_elution_profile[2 * ElutionProfileRadius + 1]; // (retention time, intensity) pairs
#endif

    QuantEntry() { memset(&pr, 0, sizeof(PrecursorEntry)); memset(fr, 0, TopF * sizeof(Fragment)); }

    template <class F> void write(F &out) { out.write((char*)this, sizeof(QuantEntry)); }

    template <class F> void read(F &in) { in.read((char*)this, sizeof(QuantEntry)); }

    friend inline bool operator < (const QuantEntry &left, const QuantEntry &right) { return left.index < right.index; }
};

class DecoyEntry {
public:
    int index = -1;
    float qvalue = 1.0, lib_qvalue = 0.0, dratio = 1.0, cscore = 0.0; // constructor below must have these intialised explicitly

    DecoyEntry() {  }
    DecoyEntry(int _index, float _qvalue, float _lib_qvalue, float _dratio, float _cscore) { index = _index; qvalue = _qvalue; lib_qvalue = _lib_qvalue; dratio = _dratio; cscore = _cscore; }

    template <class F> void write(F &out) {
        out.write((char*)this, sizeof(DecoyEntry));
    }

    template <class F> void read(F &in) {
        in.read((char*)this, sizeof(DecoyEntry));
    }

    friend inline bool operator < (const DecoyEntry &left, const DecoyEntry &right) { return left.index < right.index; }
};

class Quant {
public:
    RunStats RS;
    std::vector<QuantEntry> entries;
    std::vector<DecoyEntry> decoys; // used only for cross-run q-value calculation - for spectral library generation from DIA data;
    std::vector<int> proteins; // protein ids at <= ProteinIDQvalue
    std::vector<int> detectable_precursors; // precursors which are being searched in the particular run - bit flags in int
    double weights[pN], guide_weights[pN];
    double tandem_min, tandem_max;

    int run_index, lib_size;
    double MassAccuracy = GlobalMassAccuracy, MassAccuracyMs1 = GlobalMassAccuracyMs1;
    std::vector<double> MassCorrection, MassCorrectionMs1, MassCalSplit, MassCalSplitMs1, MassCalCenter, MassCalCenterMs1;
    double Q1Correction[2];

    template <class F> void write(F &out) {
        out.write((char*)&run_index, sizeof(int));
        out.write((char*)&lib_size, sizeof(int));
        out.write((char*)(&RS), sizeof(RunStats));
        out.write((char*)weights, pN * sizeof(double));
        out.write((char*)guide_weights, pN * sizeof(double));
        out.write((char*)&tandem_min, sizeof(double));
        out.write((char*)&tandem_max, sizeof(double));
        out.write((char*)&MassAccuracy, sizeof(double));
        out.write((char*)&MassAccuracyMs1, sizeof(double));
        out.write((char*)Q1Correction, 2 * sizeof(double));
        write_vector(out, MassCorrection);
        write_vector(out, MassCorrectionMs1);
        write_vector(out, MassCalSplit);
        write_vector(out, MassCalSplitMs1);
        write_vector(out, MassCalCenter);
        write_vector(out, MassCalCenterMs1);

        write_vector(out, proteins);
        write_vector(out, detectable_precursors);

        int size = entries.size();
        out.write((char*)&size, sizeof(int));
        for (int i = 0; i < size; i++) entries[i].write(out);

        size = decoys.size();
        out.write((char*)&size, sizeof(int));
        for (int i = 0; i < size; i++) decoys[i].write(out);

    }

    template <class F> void read_meta(F &in, int _lib_size = 0) {
        if (in.fail()) {
            qDebug() << "ERROR: cannot read the quant file\n";
            exit(0);
        }
        in.read((char*)&run_index, sizeof(int));
        in.read((char*)&lib_size, sizeof(int));
        if (_lib_size > 0 && lib_size != _lib_size) {
            qDebug() << "ERROR: a .quant file was obtained using a different spectral library / different library-free search settings\n";
            exit(0);
        }
        in.read((char*)(&RS), sizeof(RunStats));
        in.read((char*)weights, pN * sizeof(double));
        in.read((char*)guide_weights, pN * sizeof(double));
        in.read((char*)&tandem_min, sizeof(double));
        in.read((char*)&tandem_max, sizeof(double));
        in.read((char*)&MassAccuracy, sizeof(double));
        in.read((char*)&MassAccuracyMs1, sizeof(double));
        in.read((char*)Q1Correction, 2 * sizeof(double));
        read_vector(in, MassCorrection);
        read_vector(in, MassCorrectionMs1);
        read_vector(in, MassCalSplit);
        read_vector(in, MassCalSplitMs1);
        read_vector(in, MassCalCenter);
        read_vector(in, MassCalCenterMs1);
    }

    void read_meta(const std::string &file, int _lib_size = 0) {
        std::ifstream in(file, std::ifstream::binary);
        if (in.fail() || temp_folder.size()) {
            auto pos = file.find_last_of('.');
            if (pos >= 1) in = std::ifstream(location_to_file_name(file.substr(0, pos - 1)) + std::string(".quant"), std::ifstream::binary);
        }
        read_meta(in, _lib_size);
        in.close();
    }

    template <class F> void read(F &in, int _lib_size = 0) {
        read_meta(in, _lib_size);
        read_vector(in, proteins);
        read_vector(in, detectable_precursors);
        int size; in.read((char*)&size, sizeof(int));
        entries.resize(size);
        for (int i = 0; i < size; i++) entries[i].read(in);
        in.read((char*)&size, sizeof(int));
        decoys.resize(size);
        for (int i = 0; i < size; i++) decoys[i].read(in);
    }
};

std::vector<Quant> quants;


class Peptide {
public:
    int index = 0, charge = 0, length = 0, no_cal = 0;
    float mz = 0.0, iRT = 0.0, sRT = 0.0, lib_qvalue = 0.0;
    std::vector<Product> fragments;

    void init(float _mz, float _iRT, int _charge, int _index) {
        mz = _mz;
        iRT = _iRT;
        sRT = 0.0;
        charge = _charge;
        index = _index;
    }

    inline void free() {
        std::vector<Product>().swap(fragments);
    }

    template <class F> void write(F &out) {
        out.write((char*)&index, sizeof(int));
        out.write((char*)&charge, sizeof(int));
        out.write((char*)&length, sizeof(int));

        out.write((char*)&mz, sizeof(float));
        out.write((char*)&iRT, sizeof(float));
        out.write((char*)&sRT, sizeof(float));

        write_vector(out, fragments);
    }

    template <class F> void read(F &in) {
        in.read((char*)&index, sizeof(int));
        in.read((char*)&charge, sizeof(int));
        in.read((char*)&length, sizeof(int));

        in.read((char*)&mz, sizeof(float));
        in.read((char*)&iRT, sizeof(float));
        in.read((char*)&sRT, sizeof(float));

        read_vector(in, fragments);
    }

};

class FastaEntry {
public:
    std::string stripped;
    mutable std::string ids, names, genes;
    FastaEntry(const std::string &_stripped, const std::string &_ids, const std::string &_names, const std::string &_genes) { stripped = _stripped, ids = _ids, names = _names, genes = _genes; }
    friend inline bool operator < (const FastaEntry &left, const FastaEntry &right) { return left.stripped < right.stripped; }
};

bool name_included(const std::string &names, const std::string &name) {
    int start = 0, m = names.size(), n = name.size();
    if (!m || !n) return false;
    const char *p = &(name[0]);
    while (start + n <= m) {
        if (!memcmp(&(names[start]), p, n)) {
            if (start + n >= m) return true;
            if (names[start + n] == ';') return true;
        }
        while (start < m - n) {
            start++;
            if (names[start] == ';') break;
        }
        start++;
    }
    return false;
}

class Fasta {
public:
    std::string name;
    std::vector<std::string> peptides, sequences, ids;
    std::vector<std::pair<std::string, std::vector<int> > > dict;

    std::vector<FastaEntry> entries;
    std::vector<Isoform> proteins;

    void annotation_sgd(std::string &header, std::string &id, std::string &name, std::string &gene, std::string &description, bool * swissprot) {
        *swissprot = (header.find("Verified ORF") != std::string::npos);
        int start = 0, end2;
        auto end = header.find_first_of(' ');
        if (end == std::string::npos) id = fast_trim(header.substr(start + 1)), name = description = "";
        else {
            id = fast_trim(header.substr(start + 1, end - start - 1));
            for (end2 = end + 1; end2 < header.size() && header[end2] != ' '; end2++);
            end2 = header.find_first_of(' ', end2);
            name = fast_trim(header.substr(end + 1, end2 - end - 1));
            end = header.find_first_of('"', end2);
            if (end == std::string::npos) description = "";
            else {
                end2 = header.find_first_of('"', end + 1);
                if (end2 == std::string::npos) description = "";
                else description = header.substr(end + 1, end2 - end - 1);
            }
        }
        gene = name;
    }

    void annotation(std::string &header, std::string &id, std::string &name, std::string &gene, std::string &description, bool * swissprot) {
        if (header.size() < 3) {
            id = "", name = "", gene = "", description = "", *swissprot = false;
            return;
        }
        auto start = header.find("|");
        if (start == std::string::npos || (start != 3 && !(header[1] == 's' && header[2] == 'p') && !(header[1] == 't' && header[2] == 'r'))) {
            annotation_sgd(header, id, name, gene, description, swissprot);
            return;
        }
        auto end = header.find("|", start + 1);
        if (end == std::string::npos) end = header.find(" ", start + 1), name = description = "";
        else {
            auto end2 = header.find(" ", end + 1);
            if (end2 != std::string::npos) {
                name = fast_trim(header.substr(end + 1, end2 - end - 1));
                auto end3 = header.find("=", end2 + 1);
                if (end3 != std::string::npos) {
                    while (header[end3] != ' ' && end3 > end2 + 2) end3--;
                    if (end3 > end2 + 2) description = fast_trim(header.substr(end2 + 1, end3 - end2 - 1));
                    else description = "";
                } else description = fast_trim(header.substr(end2 + 1));
            } else name = fast_trim(header.substr(end + 1)), description = "";
        }
        if (end != std::string::npos) id = fast_trim(header.substr(start + 1, end - start - 1));
        else id = fast_trim(header.substr(start + 1));
        start = header.find("GN=");
        if (start == std::string::npos) gene = "";
        else {
            end = header.find(" ", start + 3);
            if (end != std::string::npos) gene = fast_trim(header.substr(start + 3, end - start - 3));
            else gene = fast_trim(header.substr(start + 3));
        }
        *swissprot = (header[1] == 's' && header[2] == 'p');
        if (SpeciesGenes && gene.size() && name.size()) {
            int pos = name.find_last_of('_');
            if (pos != std::string::npos) gene += name.substr(pos);
        }
    }

    bool load_proteins(std::vector<std::string> &files) {
        std::set<std::string> protein_ids;
        for (auto file : files) {
            if (Verbose >= 1) Time(), qDebug() << "Loading protein annotations from FASTA " << "\n";
            std::ifstream in(file);
            if (in.fail()) {
                qDebug() << "cannot read the file\n";
                return false;
            }
            bool swissprot;
            std::string line, id, name, gene, description;
            while (getline(in, line)) if (line[0] == '>') {
                    annotation(line, id, name, gene, description, &swissprot);
                    auto ins = protein_ids.insert(id);
                    if (ins.second) proteins.push_back(Isoform(id, name, gene, description, swissprot));
                }
            in.close();
        }
        std::sort(proteins.begin(), proteins.end());
        return true;
    }

    bool load(std::vector<std::string> &files, std::vector<std::string> &filter_files) {
        std::set<FastaEntry> unique;
        std::set<std::string> protein_ids, filter_peptides;
        std::vector<std::string> filters;
        std::vector<int> cuts, vec(1);
        std::map<std::string, std::vector<int> > dic;
        std::string line, sequence, peptide, id, name, gene, description, next_id, next_name, next_gene;
        bool swissprot;

        for (auto &filter : filter_files) {
            std::ifstream in(filter);
            if (!in.fail()) {
                while (getline(in, line)) if (line[0] != '>') {
                        auto trimmed = fast_trim(line);
                        trimmed.erase(std::remove(trimmed.begin(), trimmed.end(), '\t'), trimmed.end());
                        filter_peptides.insert(trimmed);
                    }
                in.close();
                if (Verbose >= 1) Time(), qDebug() << "Loaded FASTA filter " ;
            } else qDebug() << "WARNING: couldn't open the FASTA filter file " ;
            filters.insert(filters.begin(), filter_peptides.begin(), filter_peptides.end());
            filter_peptides.clear();
        }

        for (auto file : files) {

            std::ifstream in(file);
            if (in.fail()) {
                qDebug() << "cannot read the file\n";
                return false;
            }

            sequence.clear(), id.clear(), name.clear(), gene.clear(), peptide.clear();
            bool skip = false, eof = !getline(in, line);
            int skipped = 0;
            if (!eof) while (true) {
                    bool new_prot = (line.size() != 0);
                    if (new_prot) new_prot &= (line[0] == '>');
                    if (new_prot || eof) {
                        if (new_prot) {
                            annotation(line, next_id, next_name, next_gene, description, &swissprot);
                            auto ins = protein_ids.insert(next_id);
                            skip = !ins.second; // dublicate protein?
                            if (!skip) proteins.push_back(Isoform(next_id, next_name, next_gene, description, swissprot));
                            if (FastaProtDuplicates) skip = false;
                            if (skip) skipped++;
                            if (ForceSwissProt) skip = skip || !swissprot;
                        }

                        cuts.clear();
                        cuts.push_back(-1);
                        vec[0] = sequences.size();
                        for (int i = 0; i < sequence.length(); i++) {
                            char aa = sequence[i];
                            if (aa >= 'A' && aa <= 'Z') {
                                bool cut = (i == sequence.length() - 1);
                                if (!cut) cut = do_cut(aa, sequence[i + 1]);
                                if (cut) cuts.push_back(i);
                            }
                            if (i < sequence.length() - 4) {
                                auto ins = dic.insert(std::pair<std::string, std::vector<int> >(sequence.substr(i, 5), vec));
                                if (!ins.second) ins.first->second.push_back(vec[0]);
                            }
                        }
                        sequences.push_back(sequence);
                        ids.push_back(id);

                        int max_miss = MissedCleavages;
                        if (max_miss >= cuts.size()) max_miss = cuts.size() - 1;
                        for (int miss = 0; miss <= max_miss; miss++) {
                            for (int start = 0; start < cuts.size() - miss - 1; start++) {
                                bool first = (start == 0);
                                int len = cuts[start + miss + 1] - cuts[start];
                                if (len < MinPeptideLength || len > MaxPeptideLength) continue;
                                auto stripped = sequence.substr(cuts[start] + 1, cuts[start + miss + 1] - cuts[start]);
                                bool skip_peptide = false;
                                for (auto s : stripped) if (s < 'A' || s > 'Z') skip_peptide = true;
                                if (skip_peptide) continue; // handle undefined amino acids in SGD FASTA
                                if (filters.size()) {
                                    auto pos = std::lower_bound(filters.begin(), filters.end(), stripped);
                                    if (pos == filters.end()) continue;
                                    if (*pos != stripped) continue;
                                }
                                insert:
                                auto ins = unique.insert(FastaEntry(stripped, id, name, gene));
                                if (!ins.second) {
                                    if (!name_included(ins.first->ids, id)) ins.first->ids += std::string(";") + id;
                                    if (!name_included(ins.first->names, name)) ins.first->names += std::string(";") ;
                                    if (!name_included(ins.first->genes, gene)) ins.first->genes += std::string(";") + gene;
                                }

                                if (first && NMetExcision && (stripped[0] == 'M' || stripped[0] == 'X') && stripped.length() >= MinPeptideLength + 1) {
                                    stripped = stripped.substr(1);
                                    first = false;
                                    goto insert;
                                }
                            }
                        }
                        sequence.clear();
                        if (eof) break;
                        id = next_id, name = next_name, gene = next_gene;
                    } else if (!skip) sequence += line;
                    if (!getline(in, line)) eof = true, line.clear();
                }
            in.close();
            if (skipped) qDebug() << "WARNING: " << skipped << " sequences skipped due to duplicate protein ids; use --duplicate-proteins to disable skipping duplicates\n";
        }
        dict.insert(dict.begin(), dic.begin(), dic.end()); dic.clear();
        entries.insert(entries.begin(), unique.begin(), unique.end());

        std::vector<int> ind(MaxVarMods), type(MaxVarMods), mod, mod_cnt;
        for (auto &seq : unique) {
            peptide.clear();
            auto &s = seq.stripped;
            if (!FixedMods.size()) peptide = s;
            else {
                char lc = to_lower(s[0]);
                for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == lc) { // N-terminal modification: encoded with a lower-case letter for the amino acid
                            peptide += '(' + fmod.name + ')';
                            break;
                        }
                for (auto &aa : s) {
                    peptide.push_back(aa);
                    for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == aa) {
                                peptide += '(' + fmod.name + ')';
                                break;
                            }
                }
            }
            peptides.push_back(peptide);

            if (VarMods.size()) {
                int i, j, k, cnt, l = s.size(), m = Min(MaxVarMods, l);
                mod.resize(l), mod_cnt.resize(l); for (i = 0; i < l; i++) mod[i] = mod_cnt[i] = 0;
                char lc = to_lower(s[0]);
                for (int y = 0; y < VarMods.size(); y++) { // N-terminal modification: encoded with a lower-case letter for the amino acid
                    auto &vmod = VarMods[y];
                    for (int x = 0; x < vmod.aas.size(); x++) if (vmod.aas[x] == lc) {
                            mod[0] |= (1 << y), mod_cnt[0]++;
                            break;
                        }
                }
                for (i = 0; i < l; i++) {
                    for (int y = 0; y < VarMods.size(); y++) {
                        auto &vmod = VarMods[y];
                        for (int x = 0; x < vmod.aas.size(); x++) if (vmod.aas[x] == s[i]) {
                                mod[i] |= (1 << y), mod_cnt[i]++;
                                break;
                            }
                    }
                }
                for (cnt = 1; cnt <= m; cnt++) {
                    for (i = 0; i < cnt - 1; i++) ind[i] = i;
                    ind[cnt - 1] = cnt - 2;
                    while (ind[0] <= l - cnt) {
                        for (i = cnt - 1; i >= 0; i--) {
                            if (ind[i] + cnt - i < l) {
                                ind[i]++;
                                if (i < cnt - 1) for (j = i + 1; j < cnt; j++) ind[j] = ind[i] + j - i;
                                break;
                            }
                        }
                        if (i < 0) break;
                        for (i = 0; i < cnt; i++) if (!mod[ind[i]]) break;
                        if (i < cnt) continue;
                        for (i = 0; i < cnt - 1; i++) type[i] = 0;
                        type[cnt - 1] = -1;
                        while (type[0] < mod_cnt[ind[0]]) {
                            for (i = cnt - 1; i >= 0; i--) {
                                if (type[i] + 1 < mod_cnt[ind[i]]) {
                                    type[i]++;
                                    for (j = i + 1; j < cnt; j++) type[j] = 0;
                                    break;
                                }
                            }
                            if (i < 0) break;
                            peptide.clear();
                            char lc = to_lower(s[0]);
                            for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == lc) { // N-terminal modification: encoded with a lower-case letter for the amino acid
                                        peptide += '(' + fmod.name + ')';
                                        break;
                                    }
                            for (i = k = 0; i < l; i++) {
                                char aa = s[i];
                                peptide.push_back(aa);
                                if (FixedMods.size()) for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == aa) {
                                                peptide += '(' + fmod.name + ')';
                                                break;
                                            }
                                if (i > ind[cnt - 1]) continue;
                                while (k < cnt - 1 && i > ind[k]) k++;
                                if (i < ind[k]) continue;
                                int code = mod[i];
                                for (j = 0; j < type[k]; j++) code &= code - 1;
                                for (j = 0; j < 2; j++) if (code & (1 << j)) break;
                                assert(j < 2);
                                peptide += '(' + VarMods[j].name + ')';
                            }
                            peptides.push_back(peptide);
                        }
                    }
                }
            }
        }
        std::sort(proteins.begin(), proteins.end());
        return true;
    }

    void free() {
        std::vector<std::string>().swap(peptides);
        std::vector<FastaEntry>().swap(entries);
        std::vector<Isoform>().swap(proteins);
    }
};

std::vector<float> norm_totals, norm_shares, norm_ratios, norm_saved_ratios;
struct NormInfo {
    int index = 0;
    float signal = 0.0, par = 0.0;

    NormInfo() {}
    NormInfo(float _par, int _index, float _signal) { index = _index; signal = _signal; par = _par; }
    inline friend bool operator < (const NormInfo &left, const NormInfo &right) { return left.par < right.par; }
};
std::vector<NormInfo> norm_ind;

struct XIC { // for visualisation
    float RT = 0.0; // apex RT
    float qvalue = 0.0;
    float pr_mz = 0.0, fr_mz = 0.0; // masses used for extraction (i.e. mass corrected), not reference masses
    int run = -1, pr = -1, level = 0; // run index, precursor index in the library, MS level
    Product fr; // fragment info, if MS2 XIC
    std::vector<std::pair<float, float> > peaks; // RT, intensity pairs

    XIC() {}
    inline friend bool operator < (const XIC& left, const XIC& right) {
        if (left.pr < right.pr) return true;
        if (left.pr == right.pr) {
            if (left.run < right.run) return true;
            if (left.run == right.run && left.fr < right.fr) return true;
        }
        return false;
    }
    inline friend bool operator == (const XIC& left, const XIC& right) {
        return (left.pr == right.pr && left.run == right.run && !(left.fr < right.fr) && !(right.fr < left.fr));
    }

    inline int size_in_bytes() { return ((char*)&peaks) - ((char*)&RT) + sizeof(int) + peaks.size() * sizeof(std::pair<float, float>); }
};


const int fFromFasta = 1 << 0;
const int fPredictedSpectrum = 1 << 1;
const int fPredictedRT = 1 << 2;


class KARNNLIB_EXPORTS Library {
public:
    std::string name, fasta_names;
    std::vector<Isoform> proteins;
    std::vector<PG> protein_ids;
    std::vector<PG> protein_groups;
    std::vector<PG> gene_groups;
    std::vector<int> gg_index;
    std::vector<std::string> precursors; // precursor IDs in canonical format; library-based entries only
    std::vector<std::string> names;
    std::vector<std::string> genes;

    int skipped = 0;
    double iRT_min = 0.0, iRT_max = 0.0;
    bool gen_decoys = true, gen_charges = true, infer_proteotypicity = true, from_speclib = false;

    std::map<std::string, int> eg;
    std::vector<int> elution_groups; // same length as entries, indices of the elution groups
    std::vector<int> co_elution;
    std::vector<std::pair<int, int> > co_elution_index;

    class Entry {
    public:
        Lock lock;
        Library * lib;
        Peptide target;
        Peptide decoy;
        int entry_flags = 0, proteotypic = 0;
        std::string name; // precursor id
        std::set<PG>::iterator prot;
        int pid_index = 0, pg_index = 0, eg_id, best_run = -1, peak = 0, apex = 0, window = 0;
        float qvalue = 0.0, pg_qvalue = 0.0, best_fr_mz = 0.0;

        void init() {
            std::sort(target.fragments.begin(), target.fragments.end(), [](const Product &left, const Product &right) { return left.height > right.height; });
            std::sort(decoy.fragments.begin(), decoy.fragments.end(), [](const Product &left, const Product &right) { return left.height > right.height; });
            if (target.fragments.size() > MaxF) target.fragments.resize(MaxF);
            if (decoy.fragments.size() > MaxF) decoy.fragments.resize(MaxF);

            float max = 0.0;
            for (int i = 0; i < target.fragments.size(); i++) if (target.fragments[i].height > max) max = target.fragments[i].height;
            for (int i = 0; i < target.fragments.size(); i++) target.fragments[i].height /= max;

            max = 0.0;
            for (int i = 0; i < decoy.fragments.size(); i++) if (decoy.fragments[i].height > max) max = decoy.fragments[i].height;
            for (int i = 0; i < decoy.fragments.size(); i++) decoy.fragments[i].height /= max;

            target.length = decoy.length = peptide_length(name);
        }

        friend bool operator < (const Entry &left, const Entry &right) { return left.name < right.name; }

        inline void generate_decoy() {
            int i;

            decoy = target;
            auto aas = get_aas(name);
            int m = aas.size() - 2;
            if (m < 0) return;

            bool recognise = false;
            float gen_acc = 0.0;
            std::vector<Ion> pattern;
            for (auto &f : decoy.fragments) if (!(f.type & 3) || f.charge <= 0) {
                    recognise = true;
                    break;
                }
            if (recognise || ReverseDecoys) {
                pattern = recognise_fragments(name, target.fragments, gen_acc, target.charge, false, MaxRecCharge, MaxRecLoss);
                for (i = 0; i < pattern.size(); i++) {
                    if (lib->gen_charges) target.fragments[i].charge = pattern[i].charge;
                    target.fragments[i].type = (target.fragments[i].type & ~3) | pattern[i].type;
                    target.fragments[i].loss = pattern[i].loss;
                    target.fragments[i].index = pattern[i].index;
                }
                decoy.fragments = target.fragments;
                if (gen_acc > GeneratorAccuracy + E && ForceFragRec) {
                    decoy.fragments.clear();
                    lib->skipped++;
                    return;
                }
            }

            if (!ReverseDecoys) {
                int pos;
                float N_shift, C_shift;
                auto mod_state = get_mod_state(name, aas.length());

                if (!mod_state[1]) pos = 1;
                else if (!mod_state[Min(2, m + 1)]) pos = Min(2, m + 1);
                else if (!mod_state[0]) pos = 0;
                else pos = 1;
                N_shift = AA[MutateAAto[AA_index[aas[pos]]]] - AA[aas[pos]];

                if (!mod_state[m]) pos = m;
                else if (!mod_state[Max(0, m - 1)]) pos = Max(0, m - 1);
                else if (!mod_state[m + 1]) pos = m + 1;
                else pos = m;
                C_shift = AA[MutateAAto[AA_index[aas[pos]]]] - AA[aas[pos]];

                for (i = 0; i < decoy.fragments.size(); i++) {
                    double c = (double)Max(decoy.fragments[i].charge, 1);
                    if (decoy.fragments[i].type & fTypeY) decoy.fragments[i].mz += C_shift / c;
                    else decoy.fragments[i].mz += N_shift / c;
                }
            } else {
                auto seq = get_sequence(name);
                auto inv_seq = seq;

                for (i = 1; i < seq.size() - 1; i++) inv_seq[i] = seq[seq.size() - i - 1];
                for (i = 0; i < seq.size(); i++) if (Abs(inv_seq[i] - seq[i]) > 1200.0 * GeneratorAccuracy) break;
                if (i == seq.size()) inv_seq[1] += 12.0, inv_seq[seq.size() - 2] -= 12.0;

                auto dfr = generate_fragments(inv_seq, pattern, target.charge);
                for (i = 0; i < dfr.size(); i++) decoy.fragments[i].mz = dfr[i].mz;
            }
        }

        inline void annotate_charges() {
            int i;

            auto seq = get_sequence(name, &target.no_cal);
            if (!seq.size()) return;

            float gen_acc = 0.0;
            auto pattern = recognise_fragments(seq, target.fragments, gen_acc, target.charge, false, MaxRecCharge, MaxRecLoss);
            for (i = 0; i < pattern.size(); i++) {
                target.fragments[i].charge = pattern[i].charge;
                if (i < decoy.fragments.size()) decoy.fragments[i].charge = pattern[i].charge;
            }
        }

        std::vector<Ion> fragment() {
            assert(MaxF > 1000);

            int cnt = 0;
            auto seq = get_sequence(name);

            std::vector<Ion> gen;
            for (int charge = 1; charge <= 2; charge++) {
                for (int loss = loss_none; loss <= loss_NH3; loss++) {
                    auto frs = generate_fragments(seq, charge, loss, &cnt, MinFrAAs);
                    gen.insert(gen.end(), frs.begin(), frs.end());
                }
            }
            return gen;
        }

        void generate() {
            int i, cnt = 0, cnt_noloss = 0;
            auto seq = get_sequence(name);
            auto aas = get_aas(name);
            auto gen = generate_all_fragments(seq, 1, AddLosses ? loss_NH3 : loss_none, &cnt);
            std::vector<double> scores, H2O_scores, NH3_scores;
            y_scores(scores, target.charge, aas), to_exp(scores);
            if (AddLosses) y_loss_scores(H2O_scores, aas, false), y_loss_scores(NH3_scores, aas, true), to_exp(H2O_scores), to_exp(NH3_scores);
            for (i = cnt = 0; i < gen.size(); i++)
                if (gen[i].mz >= MinFrMz && gen[i].mz <= MaxFrMz) if ((gen[i].type & fTypeY) && seq.size() - gen[i].index >= MinFrAAs) cnt++;
            target.fragments.resize(cnt);
            for (i = cnt = 0; i < gen.size(); i++)
                if (gen[i].mz >= MinFrMz && gen[i].mz <= MaxFrMz) if ((gen[i].type & fTypeY) && seq.size() - gen[i].index >= MinFrAAs) {
                        target.fragments[cnt].mz = gen[i].mz;
                        if (gen[i].loss == loss_none) target.fragments[cnt].height = scores[gen[i].index];
                        else if (gen[i].loss == loss_H2O) target.fragments[cnt].height = scores[gen[i].index] * H2O_scores[gen[i].index];
                        else target.fragments[cnt].height = scores[gen[i].index] * NH3_scores[gen[i].index];
                        target.fragments[cnt].charge = gen[i].charge;
                        target.fragments[cnt].type = gen[i].type;
                        target.fragments[cnt].index = gen[i].index;
                        target.fragments[cnt].loss = gen[i].loss;
                        if (gen[i].loss == loss_none) cnt_noloss++;
                        cnt++;
                    }
            generate_decoy();
            init();
        }

        inline void free() {
            target.free();
            decoy.free();
        }

        template <class F> void write(F &out) {
            target.write(out);
            int dc = !lib->gen_decoys;
            out.write((char*)&dc, sizeof(int));
            if (dc) decoy.write(out);

            int ff = entry_flags, prt = proteotypic;
            out.write((char*)&ff, sizeof(int));
            out.write((char*)&prt, sizeof(int));
            out.write((char*)&pid_index, sizeof(int));
            write_string(out, name);
        }

        template <class F> void read(F &in) {
            target.read(in);
            int dc = 0; in.read((char*)&dc, sizeof(int));
            if (dc) decoy.read(in);

            int ff = 0, prt = 0;
            in.read((char*)&ff, sizeof(int));
            in.read((char*)&prt, sizeof(int));
            entry_flags = ff, proteotypic = prt;
            in.read((char*)&pid_index, sizeof(int));
            read_string(in, name);
        }

#if (HASH > 0)
        unsigned int hash() { return target.hash() ^ decoy.hash(); }
#endif
    };

    std::vector<Entry> entries;

    template<class F> void write(F &out) {
        int version = -1;
        out.write((char*)&version, sizeof(int));

        int gd = gen_decoys, gc = gen_charges, ip = infer_proteotypicity;
        out.write((char*)&gd, sizeof(int));
        out.write((char*)&gc, sizeof(int));
        out.write((char*)&ip, sizeof(int));

        write_string(out, name);
        write_string(out, fasta_names);
        write_array(out, proteins);
        write_array(out, protein_ids);
        write_strings(out, precursors);
        write_strings(out, names);
        write_strings(out, genes);
        out.write((char*)&iRT_min, sizeof(double));
        out.write((char*)&iRT_max, sizeof(double));
        write_array(out, entries);
        write_vector(out, elution_groups);
    }

    template<class F> void read(F &in) {
        int gd = 0, gc = 0, ip = 0, version = 0;

        in.read((char*)&version, sizeof(int));
        if (version >= 0) gd = version;
        else in.read((char*)&gd, sizeof(int));

        in.read((char*)&gc, sizeof(int));
        in.read((char*)&ip, sizeof(int));
        gen_decoys = gd, gen_charges = gc, infer_proteotypicity = ip;

        read_string(in, name);
        read_string(in, fasta_names);
        read_array(in, proteins);
        read_array(in, protein_ids);
        read_strings(in, precursors);
        read_strings(in, names);
        read_strings(in, genes);
        in.read((char*)&iRT_min, sizeof(double));
        in.read((char*)&iRT_max, sizeof(double));
        read_array(in, entries);
        for (auto &e : entries) e.lib = this;
        if (version <= -1 && in.peek() != std::char_traits<char>::eof()) read_vector(in, elution_groups);
    }

    void generate_spectra() {
        for (int i = 0; i < entries.size(); i++)
            if (entries[i].lock.set()) entries[i].generate();
    }

    void generate_all() {
        int i;
        if (Threads > 1) {
            std::vector<std::thread> threads;
            for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Library::generate_spectra, this));
            for (i = 0; i < Threads; i++) threads[i].join();
        } else generate_spectra();
        for (i = 0; i < entries.size(); i++) entries[i].lock.free();
    }

#ifdef PREDICTOR
    void smp_generate_codes(std::vector<std::pair<std::vector<long long>, int> > *_codes, std::vector<std::pair<std::string, int> > *_dict, std::vector<Lock> *_locks) {
		std::string temp_s;
		std::vector<std::string> temp_sv;
		auto &codes = *_codes;
		auto &dict = *_dict;
		auto &locks = *_locks;

		int skipped = 0;
		for (int i = 0; i < entries.size(); i++) if (locks[i].set()) {
			auto &e = entries[i];
			predictor::code_from_precursor(codes[i].first, to_charged_eg(e.name, e.target.charge), temp_s, temp_sv, dict);
			codes[i].second = i;
		}
	}

	void smp_decode_spectra(std::vector<std::vector<float> > *_spectra, std::vector<int> *_decoder, int max_frc, std::vector<Lock> *_locks) {
		auto &spectra = *_spectra;
		auto &decoder = *_decoder;
		auto &locks = *_locks;
		int i, j, cnt;

		std::vector<Ion> ions, filtered;
		for (i = 0; i < entries.size(); i++) if (locks[i].set()) {
			int ind = decoder[i];
			if (ind < 0) continue;

			auto &e = entries[i];
			auto &sp = spectra[ind];

			int tot = sp.size(), L = tot / (2 * max_frc);
			if (tot < MinGenFrNum) continue;

			ions.resize(tot); memset(&(ions[0]), 0, tot * sizeof(Ion));
			for (int c = 1; c <= max_frc; c++) for (j = 0; j < L; j++) {
				int index = L * (c - 1) + j;
				ions[index] = Ion(fTypeY, j + 1, loss_none, c, sp[index]); // trick with writing intensity in the m/z field of class Ion for later sorting
				index += L * max_frc;
				ions[index] = Ion(fTypeB, j + 3, loss_none, c, sp[index]);
			}

			std::sort(ions.begin(), ions.end(), [&](const Ion &x, const Ion &y) { return x.mz > y.mz; });
			if (ions[MinGenFrNum - 1].mz - E <= E) continue;
			double margin = Min(ions[0].mz * 0.001, ions[MinGenFrNum - 1].mz - E);
			filtered.clear(); for (auto &ion : ions) if (ion.mz > margin && (ion.type == fTypeY ? (L + 3 - ion.index) : ion.index) >= MinFrAAs) filtered.push_back(ion);
			auto frs = generate_fragments(get_sequence(e.name), filtered, e.target.charge);
			for (j = cnt = 0; j < frs.size(); j++) if (frs[j].mz >= MinFrMz && frs[j].mz <= MaxFrMz) cnt++;
			e.target.fragments.resize(Min(cnt, auxF));
			for (j = cnt = 0; j < frs.size(); j++) if (frs[j].mz >= MinFrMz && frs[j].mz <= MaxFrMz) {
				e.target.fragments[cnt++] = Product(frs[j].mz, filtered[j].mz, filtered[j].charge, filtered[j].type, filtered[j].index, loss_none);
				if (cnt >= auxF) break;
			}
			e.entry_flags |= fPredictedSpectrum;

			e.generate_decoy();
			e.init();
		}
	}

	void generate_predictions() {
		if (!entries.size()) return;
		predictor::init_predictor(1, Threads);
		P.set_instance(0);
		int i, j;

		if (Verbose >= 1) Time(), std::cout << "Encoding peptides for spectra and RTs prediction\n";
		auto dict_map = P.get_aa_indices();
		std::vector<std::pair<std::string, int> > dict(dict_map.begin(), dict_map.end());
		std::vector<std::pair<std::vector<long long>, int> > codes(entries.size());
		std::vector<Lock> locks(entries.size());

		{
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Library::smp_generate_codes, this, &codes, &dict, &locks));
			for (i = 0; i < Threads; i++) threads[i].join();
			for (auto &L : locks) L.free();
		}

		int skipped = 0;
		for (auto &c : codes) if (!c.first.size()) skipped++;
		if (skipped) cout << "WARNING: skipping " << skipped << " precursors; unrecognised modifications?\n";

		auto rt_codes = codes;
		std::sort(codes.begin(), codes.end(), [](auto &x, auto &y) { return x.first < y.first; });
		for (auto &c : rt_codes) if (c.first.size()) c.first[0] = 0;
		std::sort(rt_codes.begin(), rt_codes.end(), [](auto &x, auto &y) { return x.first < y.first; });
		std::vector<int> fr_decoder(entries.size(), -1), rt_decoder(entries.size(), -1);
		int last_fr = 0, last_rt = 0;
		auto *code = &codes[0], *rt_code = &rt_codes[0];
		for (i = 0; i < codes.size(); i++) {
			auto *curr = &codes[i];
			bool flag = (i == 0);
			if (!curr->first.size()) goto rt;
			if (!flag) if (curr->first != code->first) flag = true;
			if (!flag) curr->first.clear();
			else last_fr = i, code = curr;
			fr_decoder[curr->second] = last_fr;

		rt:
			curr = &rt_codes[i];
			if (!curr->first.size()) continue;
			flag = (i == 0);
			if (!flag) if (curr->first != rt_code->first) flag = true;
			if (!flag) curr->first.clear();
			else last_rt = i, rt_code = curr;
			rt_decoder[curr->second] = last_rt;
		}

		std::vector<std::vector<float> > spectra, rts;
		try {
			if (Verbose >= 1) Time(), std::cout << "Predicting spectra\n";
			P.predict(spectra, codes, 1, Verbose >= 3);
			if (Verbose >= 1) Time(), std::cout << "Predicting RTs\n";
			P.predict(rts, rt_codes, 2, Verbose >= 3);
		} catch (std::exception &e) { std::cout << "ERROR: " << e.what() << '\n'; }

		if (Verbose >= 1) Time(), std::cout << "Decoding predicted spectra\n";
		{
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Library::smp_decode_spectra, this, &spectra, &fr_decoder, 3, &locks));
			for (i = 0; i < Threads; i++) threads[i].join();
		}

		if (Verbose >= 1) Time(), std::cout << "Decoding RTs\n";
		for (i = 0; i < entries.size(); i++) {
			int ind = rt_decoder[i];
			if (ind >= 0) if (rts[ind].size()) {
				entries[i].target.iRT = entries[i].target.sRT = rts[ind][0] * 243.19 - 58.4;
				entries[i].entry_flags |= fPredictedRT;
			}
		}
	}
#endif

    class Info {
    public:
        Library * lib;
        int n_s; // number of samples (runs)
        int n_entries; // total number of entries
        int n_decoys; // total number of decoy entries
        int n_det_pr_block; // block size for the detectable_precursors array
        std::map<int, std::vector<std::pair<int, QuantEntry> > > map;
        std::map<int, std::vector<std::pair<int, DecoyEntry> > > decoy_map;
        std::vector<std::vector<int> > proteins;
        std::vector<int> detectable_precursors;
        std::vector<RunStats> RS;
        std::vector<std::vector<QuantEntry*> > runs;

        Info() { }

        void clear() {
            map.clear();
            decoy_map.clear();
            std::vector<std::vector<int> >().swap(proteins);
            std::vector<RunStats>().swap(RS);
            std::vector<std::vector<QuantEntry*> >().swap(runs);
        }

        void index() {
            runs.clear(); runs.resize(Max(n_s, ms_files.size()));

            for (auto it = map.begin(); it != map.end(); it++) {
                auto v = &(it->second);

                for (auto jt = v->begin(); jt != v->end(); jt++) {
                    int s = jt->first;
                    auto qe = &(jt->second);
                    runs[s].push_back(qe);
                }
            }
        }

        void load(Library * parent, std::vector<std::string> &files, std::vector<Quant> * _quants = NULL) { // files currently must be the same as ms_files
            if (Verbose >= 1) Time(), qDebug() << "Reading quantification information: " << (n_s = files.size()) << " files\n";
            clear();

            lib = parent;
            int lib_n = lib->entries.size();
            n_det_pr_block = (lib_n / sizeof(int)) + 64;
            n_entries = n_decoys = 0;
            proteins.clear();
            proteins.resize(n_s);
            detectable_precursors.resize(n_s * n_det_pr_block, 0);
            RS.resize(n_s);
            for (int i = 0; i < n_s; i++) {
                Quant Qt;
                auto Q = &Qt;
                if (_quants == NULL || !QuantInMem) {
                    std::ifstream in(std::string(files[i]) + std::string(".quant"), std::ifstream::binary);
                    if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(files[i]) + std::string(".quant"), std::ifstream::binary);
                    Qt.read(in, lib->entries.size());
                    in.close();
                } else Q = &((*_quants)[i]);
                memcpy(&(RS[i]), &(Q->RS), sizeof(RunStats));
                proteins[i] = Q->proteins;
                if (Q->lib_size != lib_n) qDebug() << "ERROR: different library used to generate the .quant file\n";
                memcpy(&detectable_precursors[i * n_det_pr_block], &Q->detectable_precursors[0], n_det_pr_block * sizeof(int));
                for (auto it = Q->entries.begin(); it != Q->entries.end(); it++) {
                    n_entries++;
                    int index = it->pr.index;
                    auto pos = map.find(index);
                    if (pos != map.end())
                        pos->second.push_back(std::pair<int, QuantEntry>(i, *it));
                    else {
                        std::vector<std::pair<int, QuantEntry> > v;
                        v.push_back(std::pair<int, QuantEntry>(i, *it));
                        map.insert(std::pair<int, std::vector<std::pair<int, QuantEntry> > >(index, v));
                    }
                }
                if (_quants == NULL || !QuantInMem) std::vector<QuantEntry>().swap(Q->entries);

                for (auto it = Q->decoys.begin(); it != Q->decoys.end(); it++) {
                    n_decoys++;
                    int index = it->index;
                    auto pos = decoy_map.find(index);
                    if (pos != decoy_map.end())
                        pos->second.push_back(std::pair<int, DecoyEntry>(i, DecoyEntry(index, it->qvalue, it->lib_qvalue, it->dratio, it->cscore)));
                    else {
                        std::vector<std::pair<int, DecoyEntry> > v;
                        v.push_back(std::pair<int, DecoyEntry>(i, DecoyEntry(index, it->qvalue, it->lib_qvalue, it->dratio, it->cscore)));
                        decoy_map.insert(std::pair<int, std::vector<std::pair<int, DecoyEntry> > >(index, v));
                    }
                }
            }

            index();
        }

        void load(Library * parent, Quant &quant) {
            clear();

            lib = parent;
            int lib_n = lib->entries.size();
            n_det_pr_block = (lib_n / sizeof(int)) + 64;
            n_entries = 0;
            proteins.clear();
            proteins.resize(n_s = 1);
            detectable_precursors.resize(n_det_pr_block, 0);

            Quant &Q = quant;
            proteins[0] = Q.proteins;
            if (Q.lib_size != lib_n) qDebug() << "ERROR: different library used to generate the .quant file\n";
            memcpy(&detectable_precursors[0], &Q.detectable_precursors[0], n_det_pr_block * sizeof(int));
            for (auto it = Q.entries.begin(); it != Q.entries.end(); it++) {
                n_entries++;
                int index = it->pr.index;
                auto pos = map.find(index);
                if (pos != map.end())
                    pos->second.push_back(std::pair<int, QuantEntry>(0, *it));
                else {
                    std::vector<std::pair<int, QuantEntry> > v;
                    v.push_back(std::pair<int, QuantEntry>(0, *it));
                    map.insert(std::pair<int, std::vector<std::pair<int, QuantEntry> > >(index, v));
                }
            }

            for (auto it = Q.decoys.begin(); it != Q.decoys.end(); it++) {
                n_decoys++;
                int index = it->index;
                auto pos = decoy_map.find(index);
                if (pos != decoy_map.end())
                    pos->second.push_back(std::pair<int, DecoyEntry>(0, DecoyEntry(index, it->qvalue, it->lib_qvalue, it->dratio, it->cscore)));
                else {
                    std::vector<std::pair<int, DecoyEntry> > v;
                    v.push_back(std::pair<int, DecoyEntry>(0, DecoyEntry(index, it->qvalue, it->lib_qvalue, it->dratio, it->cscore)));
                    decoy_map.insert(std::pair<int, std::vector<std::pair<int, DecoyEntry> > >(index, v));
                }
            }

            index();
        }

        void quantify() {
            int i, j, k, m, pos;

            if (Verbose >= 1) Time(), qDebug() << "Quantifying peptides\n";
            std::vector<std::pair<float, int> > fr_score(1000), fr_ordered(1000);
            std::vector<bool> exclude(1024);
            for (auto it = map.begin(); it != map.end(); it++) {
                auto v = &(it->second);
                int pr_index = it->first;
                auto &lib_e = lib->entries[pr_index].target;

                fr_score.clear();
                for (auto jt = (*v).begin(); jt != (*v).end(); jt++) if (jt->second.pr.qvalue <= MaxQuantQvalue) {
                        for (i = 0; i < jt->second.fr_n; i++) {
                            int ind = jt->second.fr[i].index;
                            if (ind >= fr_score.size()) fr_score.resize(ind + 1, std::pair<float, int>(0.0, 0));
                            fr_score[ind].first += jt->second.fr[i].corr;
                            fr_score[ind].second++;
                        }
                    }
                for (auto &s : fr_score) if (s.second) s.first /= (double)s.second;

                if (ExcludeSharedFragments) { // exclude (from quantification) fragments shared by peptides in the same elution group
                    int eg = lib->elution_groups[pr_index];
                    auto &ce = lib->co_elution_index[eg];
                    for (int pos = ce.first; pos < ce.first + ce.second; pos++) {
                        int next = lib->co_elution[pos];
                        if (next != pr_index && lib->entries[next].target.charge == lib->entries[pr_index].target.charge) { // difference only in the label
                            auto &lf = lib_e;
                            auto &ls = lib->entries[next].target;
                            for (i = 0; i < lf.fragments.size(); i++) {
                                double margin = GeneratorAccuracy * lf.fragments[i].mz;
                                bool flag = false;
                                for (j = 0; j < ls.fragments.size(); j++)
                                    if (Abs(lf.fragments[i].mz - ls.fragments[j].mz) < margin) {
                                        flag = true;
                                        break;
                                    }
                                if (flag) fr_score[i].first -= 2.0;
                            }
                        }
                    }
                }

                if (RestrictFragments) { // exclude fragments from quantification based on their library annotation
                    auto &frs = lib_e.fragments;
                    assert(fr_score.size() <= frs.size());
                    for (j = 0; j < fr_score.size(); j++) if (frs[j].type & fExclude) fr_score[j].first -= INF;
                }

                fr_ordered.resize(fr_score.size());
                fr_ordered.assign(fr_score.begin(), fr_score.end());
                std::sort(fr_ordered.begin(), fr_ordered.end());
                for (pos = fr_ordered.size() - 3; pos < fr_ordered.size() - 1; pos++) if (fr_ordered[pos].first > -1.0 + E) break;
                double margin = Max(-INF / 2.0, fr_ordered[pos].first - E);
                if (NoFragmentSelectionForQuant) margin = -INF / 2.0;
                if (GenFrExclusionInfo) for (j = 0; j < fr_score.size(); j++) if (fr_score[j].second) {
                            auto &fr = lib_e.fragments[j];
                            if (fr_score[j].first >= margin) fr.type &= ~fExclude;
                            else fr.type |= fExclude;
                        }

                for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                    for (i = 0, jt->second.pr.quantity = jt->second.pr.ratio = jt->second.pr.quality = 0.0; i < jt->second.fr_n; i++) {
                        int ind = jt->second.fr[i].index;
                        if (fr_score[ind].first >= margin) {
                            jt->second.pr.quantity += jt->second.fr[i].quantity[qFiltered];
                            jt->second.pr.quality += jt->second.fr[i].corr * jt->second.fr[i].quantity[qFiltered];
                        }
                    }
                    if (jt->second.pr.quantity > E) jt->second.pr.quality /= jt->second.pr.quantity;
                    jt->second.pr.norm = jt->second.pr.level = jt->second.pr.quantity;
                }
            }

            if (TranslatePeaks) { // calculate ratios between labelled and unlabelled peptides; requires TranslatePeaks to ensure co-elution
                for (auto it = map.begin(); it != map.end(); it++) {
                    auto v = &(it->second);
                    int pr_index = it->first;
                    auto &lib_e = lib->entries[pr_index].target;

                    if (lib->entries[pr_index].entry_flags & fFromFasta) continue; // only spectral library entries

                    int eg = lib->elution_groups[pr_index];
                    auto &ce = lib->co_elution_index[eg];

                    exclude.clear();
                    if (ExcludeSharedFragments) { // exclude (from quantification) fragments shared by peptides in the same elution group
                        int eg = lib->elution_groups[pr_index];
                        auto &ce = lib->co_elution_index[eg];
                        for (int pos = ce.first; pos < ce.first + ce.second; pos++) {
                            int next = lib->co_elution[pos];
                            if (next != pr_index && lib->entries[next].target.charge == lib->entries[pr_index].target.charge) { // difference only in the label
                                auto &lf = lib_e;
                                auto &ls = lib->entries[next].target;
                                for (i = 0; i < lf.fragments.size(); i++) {
                                    double margin = GeneratorAccuracy * lf.fragments[i].mz;
                                    bool flag = false;
                                    for (j = 0; j < ls.fragments.size(); j++)
                                        if (Abs(lf.fragments[i].mz - ls.fragments[j].mz) < margin) {
                                            flag = true;
                                            break;
                                        }
                                    if (flag) {
                                        if (i >= exclude.size()) exclude.resize(i + 1, false);
                                        exclude[i] = true;
                                    }
                                }
                            }
                        }
                    }

                    const int fr_n = 3;
                    std::vector<float> fratio(fr_n);
                    for (int pos = ce.first; pos < ce.first + ce.second; pos++) {
                        int next = lib->co_elution[pos];
                        if (next != pr_index) if (lib->entries[next].target.charge == lib->entries[pr_index].target.charge) {
                                auto p = map.find(next);
                                if (p != map.end()) if (p->first == next) {
                                        auto nv = &(p->second);
                                        auto kt = (*nv).begin();
                                        for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                                            for (; kt != (*nv).end(); kt++) if (kt->first >= jt->first) break;
                                            if (kt->first == jt->first) {
                                                fr_score.clear();
                                                for (i = 0; i < jt->second.fr_n; i++) {
                                                    int ind = jt->second.fr[i].index;
                                                    if (ind >= fr_score.size()) fr_score.resize(ind + 1, std::pair<float, int>(0.0, 0));
                                                    fr_score[ind].first = jt->second.fr[i].corr, fr_score[ind].second = lib->entries[pr_index].target.fragments[ind].ion_code();
                                                }
                                                for (i = 0; i < kt->second.fr_n; i++) {
                                                    int ind = kt->second.fr[i].index;
                                                    if (ind >= fr_score.size()) fr_score.resize(ind + 1, std::pair<float, int>(0.0, 0));
                                                    if (fr_score[ind].second != lib->entries[next].target.fragments[ind].ion_code()) fr_score[i].first = -INF;
                                                    else if (kt->second.fr[i].corr < fr_score[ind].first) fr_score[ind].first = kt->second.fr[i].corr;
                                                }
                                                for (i = 0; i < fr_score.size(); i++) {
                                                    fr_score[i].second = i;
                                                    if (i < exclude.size()) if (exclude[i]) fr_score[i].first -= INF * 0.5;
                                                }
                                                std::sort(fr_score.begin(), fr_score.end());
                                                if (fr_score.size()) if (fr_score[fr_score.size() - 1].first > E) {
                                                        int fcnt = 0;
                                                        float f1 = 0.0, f2 = 0.0, f1f[fr_n], f2f[fr_n];
                                                        for (i = 0; i < fr_n; i++) f1f[i] = f2f[i] = 0.0;
                                                        for (int ifr = fr_score.size() - 1; ifr >= 0 && ifr >= fr_score.size() - fr_n; ifr--) {
                                                            if (ifr < fr_score.size() - 1) if (fr_score[ifr].first < 0.8) break;
                                                            int fr = fr_score[ifr].second;
                                                            for (i = 0; i < jt->second.fr_n; i++) if (jt->second.fr[i].index == fr) {
                                                                    f1f[fcnt] += jt->second.fr[i].quantity[qFiltered];
                                                                    break;
                                                                }
                                                            for (i = 0; i < kt->second.fr_n; i++) if (kt->second.fr[i].index == fr) {
                                                                    f2f[fcnt] += kt->second.fr[i].quantity[qFiltered];
                                                                    break;
                                                                }
                                                            fcnt++;
                                                            if (fcnt >= fr_n) break;
                                                        }
                                                        fratio.resize(fr_n);
                                                        for (i = fcnt = 0; i < fr_n; i++) if (f1f[i] > E) fratio[fcnt++] = f2f[i] / f1f[i];
                                                        if (fcnt <= 2) {
                                                            for (i = 0; i < fr_n; i++) f1 += f1f[i], f2 += f2f[i];
                                                            if (f1 > E) jt->second.pr.ratio += (jt->second.pr.quantity * f2) / f1;
                                                        } else {
                                                            fratio.resize(fcnt); std::sort(fratio.begin(), fratio.end());
                                                            jt->second.pr.ratio += jt->second.pr.quantity * fratio[fcnt / 2];
                                                        }
                                                        continue;
                                                    }
                                                jt->second.pr.ratio += kt->second.pr.quantity;
                                            }
                                        }
                                    }
                            }
                    }
                    for (auto jt = (*v).begin(); jt != (*v).end(); jt++) jt->second.pr.ratio = (jt->second.pr.ratio > E ? (jt->second.pr.quantity / jt->second.pr.ratio) : -1.0);
                }
            }

            if (ms_files.size() <= 1 || IndividualReports) return;

            // simple normalisation by the total signal first
            double av = 0.0;
            std::vector<float> sums(ms_files.size());
            for (auto it = map.begin(); it != map.end(); it++) {
                auto v = &(it->second);
                for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                    int index = jt->first;
                    if (jt->second.pr.qvalue <= NormalisationQvalue) sums[index] += jt->second.pr.quantity;
                }
            }
            for (i = k = 0; i < sums.size(); i++) if (sums[i] > E) av += sums[i], k++;
            if (k) av /= (double)k;
            if (av <= E) {
                Warning("not enough peptides for normalisation");
                return;
            }

            for (i = 0; i < sums.size(); i++) if (sums[i] <= E) sums[i] = av;
            for (auto it = map.begin(); it != map.end(); it++) {
                auto v = &(it->second);
                for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                    int index = jt->first;
                    jt->second.pr.level = (jt->second.pr.quantity * av) / sums[index];
                }
            }

            // advanced normalisation
            std::vector<float> score(map.size());
            std::vector<float> x(ms_files.size());
            i = m = 0;
            for (auto it = map.begin(); it != map.end(); it++, i++) {
                auto v = &(it->second);

                k = 0;
                memset(&x[0], 0, x.size() * sizeof(float));
                for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                    int index = jt->first;
                    if (jt->second.pr.qvalue <= NormalisationQvalue) x[index] = jt->second.pr.level, k++, m++;
                }
                if (k >= 2) {
                    double u = mean(&(x[0]), x.size());
                    if (u > E) score[i] = sqrt(var(&(x[0]), x.size())) / u;
                    else score[i] = INF;
                }
                else score[i] = INF;
            }
            if (!m) {
                Warning("not enough peptides for normalisation");
                return;
            }

            auto ordered = score;
            std::sort(ordered.begin(), ordered.end());
            int average_number = m / ms_files.size();
            int used = Max(1, int(NormalisationPeptidesFraction * (double)average_number));
            double margin = ordered[used] + E;

            for (i = 0; i < sums.size(); i++) sums[i] = 0.0;
            i = 0;
            for (auto it = map.begin(); it != map.end(); it++, i++) {
                auto v = &(it->second);
                if (score[i] > margin) continue;

                for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                    int index = jt->first;
                    if (jt->second.pr.qvalue <= NormalisationQvalue) sums[index] += jt->second.pr.level;
                }
            }

            for (i = k = 0, av = 0.0; i < sums.size(); i++) if (sums[i] > E) av += sums[i], k++;
            if (k) av /= (double)k;
            if (av < E) {
                Warning("cannot perform normalisation");
                return;
            }
            for (i = 0; i < sums.size(); i++) if (sums[i] <= E) sums[i] = av;
            for (auto it = map.begin(); it != map.end(); it++) {
                auto v = &(it->second);

                for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                    int index = jt->first;
                    jt->second.pr.norm = (jt->second.pr.level * av) / sums[index];
                }
            }

            if (LocalNormalisation) for (int iter = 0; iter < 2; iter++) { // iter = 0: RT-dependent, iter = 1: signal-dependent normalisation
                    if (iter == 0 && NoRTDepNorm) continue;
                    if (iter == 1 && NoSigDepNorm) continue;

                    norm_totals.clear(), norm_shares.clear();
                    norm_totals.resize(map.size(), 0.0), norm_shares.resize(map.size(), 0.0);
                    std::vector<int> map_index(lib->entries.size(), 0);

                    i = 0;
                    for (auto it = map.begin(); it != map.end(); it++, i++) {
                        auto v = &(it->second);
                        int cnt = 0;
                        map_index[it->first] = i;
                        for (auto jt = (*v).begin(); jt != (*v).end(); jt++)
                            if (jt->second.pr.qvalue <= NormalisationQvalue) norm_totals[i] += jt->second.pr.norm, cnt++;
                        if (cnt) norm_shares[i] = 1.0 / (double)cnt;
                    }

                    norm_ind.resize(map.size()), norm_saved_ratios.resize(map.size());
                    norm_ratios.resize(2 * LocNormRadius + 1);

                    for (k = 0; k < ms_files.size(); k++) {
                        norm_ind.clear(), norm_saved_ratios.clear();

                        assert(runs.size() > k);
                        for (auto it = runs[k].begin(); it != runs[k].end(); it++) {
                            auto jt = *it;
                            i = map_index[jt->pr.index];
                            if (jt->pr.qvalue <= NormalisationQvalue && norm_shares[i] > E) if (norm_shares[i] * (double)ms_files.size() <= 2.0)
                                    norm_ind.push_back(NormInfo(iter == 0 ? jt->pr.RT : jt->pr.norm, i, jt->pr.norm));
                        }

                        std::sort(norm_ind.begin(), norm_ind.end());
                        norm_saved_ratios.resize(norm_ind.size(), 0.0);

                        for (auto it = runs[k].begin(); it != runs[k].end(); it++) {
                            auto jt = *it;

                            float ratio = 1.0;
                            auto pos = std::lower_bound(norm_ind.begin(), norm_ind.end(), NormInfo(iter == 0 ? jt->pr.RT : jt->pr.norm, 0, 0.0));
                            int ind = std::distance(norm_ind.begin(), pos), low = Max(0, ind - LocNormRadius), high = low + 2 * LocNormRadius;
                            if (high > norm_ind.size()) high = norm_ind.size(), low = Max(0, high - 2 * LocNormRadius);
                            if (norm_saved_ratios[Min(ind, norm_ind.size() - 1)] > E) ratio = norm_saved_ratios[Min(ind, norm_ind.size() - 1)];
                            else if (high - low <= Max(16, LocNormRadius / 16)) ratio = 1.0;
                            else {
                                norm_ratios.clear();
                                for (j = low; j < high; j++) {
                                    float tot = norm_totals[norm_ind[j].index] * norm_shares[norm_ind[j].index];
                                    float r = tot / Max(E, norm_ind[j].signal);
                                    norm_ratios.push_back(r);
                                }
                                std::sort(norm_ratios.begin(), norm_ratios.end());
                                if (norm_ratios.size() & 1) ratio = norm_ratios[norm_ratios.size() / 2];
                                else ratio = 0.5 * (norm_ratios[(norm_ratios.size() / 2) - 1] + norm_ratios[norm_ratios.size() / 2]);
                                norm_saved_ratios[Min(ind, norm_ind.size() - 1)] = ratio = Min(LocNormMax, Max(1.0 / LocNormMax, ratio));
                            }
                            jt->pr.norm *= ratio;
                        }
                    }
                }
        }
    };

    Info info;

    void extract_proteins() {
        std::set<Isoform> prot;
        std::string word;

        for (auto &pg : protein_ids) {
            std::stringstream list(pg.ids);
            while (std::getline(list, word, ';')) {
                auto ins = prot.insert(Isoform(fast_trim(word)));
                for (auto &pr : pg.precursors) ins.first->precursors.insert(pr);
            }
        }

        proteins.clear();
        proteins.insert(proteins.begin(), prot.begin(), prot.end());
        prot.clear();

        for (int i = 0; i < proteins.size(); i++) {
            auto &p = proteins[i];
            for (auto &pr : p.precursors) {
                auto &pid = protein_ids[entries[pr].pid_index];
                pid.proteins.insert(i);
            }
        }

        if (infer_proteotypicity) {
            if (Verbose >= 1 && (!FastaSearch || GuideLibrary)) Time(), qDebug() << "Finding proteotypic peptides (assuming that the list of UniProt ids provided for each peptide is complete)\n";
            for (auto &e : entries) if (protein_ids[e.pid_index].proteins.size() <= 1) e.proteotypic = true;
        }
    }

    void assemble_elution_groups() {
        if (Verbose >= 1) Time(), qDebug() << "Assembling elution groups\n";
        eg.clear(), elution_groups.clear(), elution_groups.reserve(entries.size());
        for (auto &e : entries) {
            auto name = to_eg(e.name);
            auto egp = eg.insert(std::pair<std::string, int>(name, eg.size()));
            elution_groups.push_back(e.eg_id = egp.first->second);
        }
        eg.clear();
    }

    void elution_group_index() {
        int i, egt = 0, peg;
        for (i = 0; i < elution_groups.size(); i++) if (elution_groups[i] > egt) egt = elution_groups[i]; egt++;
        co_elution_index.resize(egt, std::pair<int, int>(0, 1));
        std::set<std::pair<int, int> > ce;
        for (i = 0; i < elution_groups.size(); i++) ce.insert(std::pair<int, int>(elution_groups[i], i));
        co_elution.resize(ce.size()); i = peg = 0;
        if (co_elution_index.size()) co_elution_index[0].second = 0;
        for (auto it = ce.begin(); it != ce.end(); it++, i++) {
            co_elution[i] = it->second;
            if (it->first != peg) co_elution_index[it->first].first = i;
            else co_elution_index[it->first].second++;
            peg = it->first;
        }
        ce.clear();
    }

    bool load_sptxt(const char * file_name, double acc = 0.000005) {
        std::string line, aas, spectrum_name, pep_name, name, ann, mod_name, comment;
        std::map<std::string, Entry> map;
        std::set<PG> prot; prot.insert(PG(""));
        Entry e; e.lib = this;
        float pr_mz = 0.0, pr_rt = 0.0;
        bool first = false, skip = true, omit = true, no_rt = true;
        int i, j, charge = 0, no_rt_cnt = 0;
        gen_charges = false;

        qDebug() << "WARNING: support for .sptxt/.msp spectral libraries is experimental; fragment ions must be annotated; "
              << acc * 1000000.0 << " ppm mass accuracy filtering will be used; use --sptxt-acc <ppm> to change this setting\n";

        auto ins = map.insert(std::pair<std::string, Entry>("", e)).first;
        map.clear();

        std::ifstream sp(file_name, std::ifstream::in);
        while (std::getline(sp, line)) {
            if (line.size() >= 6 && (line[0] < '0' || line[0] > '9')) {
                if (!memcmp(&(line[0]), "Name:", 5)) {
                    spectrum_name = trim(line.substr(5));
                    auto cpos = spectrum_name.find_last_of('/');
                    if (cpos == std::string::npos) omit = true;
                    else {
                        pep_name = spectrum_name.substr(0, cpos);
                        charge = std::stoi(&(spectrum_name[cpos + 1]));
                        name = to_canonical(pep_name, charge);
                        omit = false;
                    }
                } else if (line.size() >= 15 && !omit) {
                    if (!memcmp(&(line[0]), "PrecursorMZ:", 12))
                        pr_mz = to_double(&(line[12])), first = true, skip = false;
                    else if (!memcmp(&(line[0]), "Comment", 7)) {
                        comment = line;
                        int ppos = line.find("Parent=");
                        if (ppos != std::string::npos) {
                            pr_mz = to_double(&(line[ppos + 7])), first = true, skip = false, pr_rt = 0.0, no_rt = true;

                            int rtpos = line.find("RT=");
                            if (rtpos != std::string::npos) pr_rt = to_double(&(line[rtpos + 3])), no_rt = false;
                            else {
                                rtpos = line.find("RetentionTime=");
                                if (rtpos != std::string::npos) pr_rt = to_double(&(line[rtpos + 14])), no_rt = false;
                            }
                            if (no_rt) no_rt_cnt++;

                            int mpos = line.find("Mods=");
                            if (mpos != std::string::npos) if (line[mpos + 5] != '0') {
                                    aas = get_aas(pep_name); name.clear();
                                    int mods = std::stoi(&(line[mpos + 5])), curr_aa = 0, nterm = 0;
                                    auto pos = line.find('/', mpos + 5) + 1;
                                    for (i = 0; i < mods; i++) {
                                        int mod_pos = std::stoi(&(line[pos]));
                                        if (mod_pos == -1) mod_pos = 0, nterm = 1;
                                        if (aas.size() <= mod_pos) {
                                            qDebug() << "ERROR: modification position " << mod_pos << " exceeds the number of amino acids in the precursor " << ", see ";
                                            return false;
                                        }
                                        pos = line.find(',', pos + 1) + 1;
                                        if (aas[mod_pos] != line[pos]) {

                                            return false;
                                        }
                                        mod_name.clear();
                                        pos = line.find(',', pos + 1) + 1;
                                        for (j = pos; line[j] && line[j] != '/' && line[j] != ' '; j++) mod_name += line[j];
                                        if (i < mods - 1) if (line[j] != '/') {

                                                return false;
                                            }
                                        pos = j + 1;

                                        if (curr_aa) nterm = 0;
                                        if (nterm) name += '(' + mod_name + ')';
                                        for (; curr_aa <= mod_pos; curr_aa++) name += aas[curr_aa];
                                        if (!nterm) name += '(' + mod_name + ')';
                                    }
                                    for (; curr_aa < aas.size(); curr_aa++) name += aas[curr_aa];
                                    name += std::to_string(charge);
                                }
                        }
                    }
                }
            }

            if (!skip) if (line.size() >= 4) {
                    if (line[0] >= '0' && line[0] <= '9') {
                        float mz = to_double(line), ref = 0.0;
                        auto pos = line.find_first_of('\t');
                        if (pos != std::string::npos && pos < line.size()) ref = to_double(&(line[++pos]));
                        if (ref > 0.0) {
                            pos = line.find('\t', pos) + 1;
                            if (line[pos] == '\"') pos++;
                            if (line[pos] == 'y' || line[pos] == 'b') {
                                ann = line.substr(pos);
                                auto apos = ann.find_first_of('/');
                                if (apos != std::string::npos) {
                                    double err = to_double(&(ann[apos + 1]));
                                    if (Abs(err / mz) < acc) {
                                        auto is_i = ann.find_first_of('i');
                                        if (is_i == std::string::npos || is_i > apos) {
                                            int type = (line[pos] == 'y' ? fTypeY : fTypeB);
                                            int num = std::stoi(&(line[pos + 1]));
                                            if (type == fTypeY) num = peptide_length(pep_name) - num;
                                            int frc = 1;
                                            auto frcpos = ann.find_first_of('^');
                                            if (frcpos != std::string::npos && frcpos < apos) frc = std::stoi(&(ann[frcpos + 1]));
                                            int loss = loss_none;
                                            auto lpos = ann.find_first_of('-');
                                            if (lpos != std::string::npos && lpos < apos) {
                                                loss = std::stoi(&(ann[lpos + 1]));
                                                if (loss == 17) loss = loss_NH3;
                                                else if (loss = 18) loss = loss_H2O;
                                                else if (loss = 28) loss = loss_CO;
                                                else loss = loss_other;
                                            }

                                            Product p(mz, ref, frc, type, num, loss);
                                            if (first) {
                                                first = false;
                                                auto inserted = map.insert(std::pair<std::string, Entry>(name, e));
                                                if (inserted.second) {
                                                    ins = inserted.first;
                                                    ins->second.name = ins->first;
                                                    ins->second.prot = prot.begin();
                                                } else {
                                                    skip = true;

                                                }
                                            }

                                            auto pep = &(ins->second.target);
                                            pep->mz = pr_mz;
                                            pep->iRT = no_rt ? 0.0 : pr_rt;
                                            pep->charge = charge;
                                            pep->fragments.push_back(p);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
        }
        sp.close();
        eg.clear();

        if (no_rt_cnt) qDebug() << "WARNING: " << no_rt_cnt << " precursors without retention time information\n";

        entries.resize(map.size());
        precursors.resize(map.size());
        i = 0;
        for (auto it = map.begin(); it != map.end(); it++, i++) {
            precursors[i] = it->first;
            entries[i] = it->second;
            entries[i].target.index = entries[i].decoy.index = i;
            if (entries[i].target.iRT < iRT_min) iRT_min = entries[i].target.iRT;
            if (entries[i].target.iRT > iRT_max) iRT_max = entries[i].target.iRT;
            it->second.prot->precursors.push_back(i);
        }
        protein_ids.insert(protein_ids.begin(), prot.begin(), prot.end());

        return true;
    }

    bool load(const char * file_name) {
        int colInd[libCols];

        name = std::string(file_name);
        if (Verbose >= 1) Time(), qDebug() << "Loading spectral library " << "\n";

        if (get_extension(name) == std::string(".speclib")) {

            std::ifstream speclib(file_name, std::ifstream::binary);

            if (speclib.fail()) {
                qDebug() << "cannot read the file\n";
                return false;
            }
            read(speclib);
            speclib.close();

            from_speclib = true;
            if (PGLevel != PGLevelSet) {
                if (PGLevelSet == 2 && genes.size() >= 2) PGLevel = 2;
                else if (PGLevelSet == 1 && names.size() >= 2) PGLevel = 1;
            }
            if (fasta_names.size()) if (Verbose >= 1) Time(), qDebug() << "Library annotated with sequence database(s): "  << "\n";
            if (!fasta_files.size() && (genes.size() >= 2 || names.size() >= 2)) library_protein_stats();
            if (!elution_groups.size()) assemble_elution_groups();

            int gc = gen_charges;
            gen_charges = false;
            for (auto &e : entries) {
                e.entry_flags &= ~fFromFasta;
                if (gc) for (auto &f : e.target.fragments) if (!(int)f.charge) gen_charges = true, gc = 0;
            }
        }

        else if (get_extension(name) == std::string(".sptxt") || get_extension(name) == std::string(".msp")) {
            load_sptxt(file_name, SptxtMassAccuracy);
            assemble_elution_groups();
        }

        else {
            std::ifstream csv(file_name, std::ifstream::in);
            if (csv.fail()) {
                qDebug() << "cannot read the file\n";
                return false;
            }

            int i, cnt = 0;
            bool ftw = false, fcw = false, ftmw = false, flw = false, fragment_num_info = false, fragment_loss_info = false, fragment_type = false;
            std::map<std::string, Entry> map;
            Entry e; e.lib = this;

            std::string line, prev_id = "", id;
            std::vector<std::string> words(1000);
            char delim = '\t';
            if (std::string(file_name).find(".csv") != std::string::npos) delim = ','; // not for OpenSWATH libraries
            auto ins = map.insert(std::pair<std::string, Entry>("", e)).first; // insert a dummy entry to create a variable (ins) of the desired type
            map.clear(); // remove the dummy entry

            std::set<PG> prot;
            eg.clear();
            while (std::getline(csv, line)) {
                if (!cnt) {
                    if (delim == '\t' && line.find("\t") == std::string::npos) delim = ',';
                    if (delim == ',' && line.find(",") == std::string::npos) delim = '\t';
                }
                int in_quote = 0, nw = 0, start = 0, end = 0;
                for (; end < line.size(); end++) {
                    if (line[end] == delim && !in_quote) {
                        int en = end, st = start;
                        if (line[st] == '\"') st++, en--;
                        if (en > st) {
                            if (nw > words.size()) words.resize(words.size() * 2 + 1);
                            words[nw].resize(en - st);
                            for (i = st; i < en; i++) words[nw][i - st] = line[i];
                        } else words[nw].clear();
                        start = end + 1, nw++;
                    } else if (line[end] == '\"') in_quote ^= 1;
                }
                if (line[start] == '\"') start++, end--;
                if (end > start) {
                    if (nw > words.size()) words.resize(words.size() * 2 + 1);
                    words[nw].resize(end - start);
                    for (i = start; i < end; i++) words[nw][i - start] = line[i];
                    nw++;
                }
                cnt++;

                if (cnt == 1) { // header
                    std::vector<std::string>::iterator it, loc;
                    for (i = 0; i < libCols; i++) colInd[i] = -1;
                    for (i = 0; i < nw; i++) {
                        for (auto jt = library_headers.begin(); jt != library_headers.end(); jt++)
                            if (jt->find(words[i]) != std::string::npos) {
                                int ind = std::distance(library_headers.begin(), jt);
                                if (colInd[ind] < 0) colInd[ind] = i;
                                else if (words[colInd[ind]].size() < words[i].size()) colInd[ind] = i;
                                break;
                            }
                    }
                    for (i = 0; i < libCols; i++) if (colInd[i] < 0 && i < libPID) {
                            if (i == 0) {
                                for (int j = 0; j < nw; j++) if (trim(words[j]) == std::string("ModifiedPeptideSequence")) {
                                        colInd[0] = j;
                                        break;
                                    }
                                if (colInd[0] >= 0) continue;
                            }
                            if (Verbose >= 1) qDebug() << "WARNING: cannot find column " + QString::fromStdString(library_headers[i]) << "\n";
                            csv.close();
                            return false;
                        }
                    if (colInd[libPID] < 0)
                        for (i = 0; i < nw; i++)
                            if (words[i] == std::string("ProteinId") || words[i] == std::string("ProteinID") || words[i] == std::string("Protein")) {
                                colInd[libPID] = i;
                                break;
                            }

                    if (colInd[libPT] >= 0) infer_proteotypicity = false;
                    if (colInd[libFrCharge] >= 0) {
                        gen_charges = false;
                        if (colInd[libFrType] >= 0) {
                            fragment_type = true;
                            if (colInd[libFrNumber] >= 0) {
                                fragment_num_info = true;
                                if (colInd[libFrLoss] >= 0) fragment_loss_info = true;
                                else qDebug() << "WARNING: no neutral loss information found in the library - assuming fragments without losses\n";
                            } else qDebug() << "WARNING: no fragment number information found in the library\n";
                        } else qDebug() << "WARNING: no fragment type information found in the library\n";
                    } else qDebug() << "WARNING: no fragment charge information found in the library - assuming fragments with charge 1\n";
                    continue;
                }

                Product p(to_double(words[colInd[libFrMz]]), to_double(words[colInd[libFrI]]), (!gen_charges ? std::stoi(words[colInd[libFrCharge]]) : 1));
                if (p.charge <= 0) {
                    p.charge = 1;
                    if (!fcw) fcw = true, std::cout << "WARNING: fragment charge 0 specified - assuming 1\n";
                }
                if (fragment_type) {
                    auto &wt = words[colInd[libFrType]];
                    if (wt.size()) {
                        if (wt[0] == 'y') p.type = fTypeY;
                        else if (wt[0] == 'b') p.type = fTypeB;

                    } else if (!ftmw) ftmw = true, p.type = 0, qDebug() << "WARNING: fragment type missing for row number " << cnt << "\n";
                    if (p.type && fragment_num_info) {
                        p.index = std::stoi(words[colInd[libFrNumber]]);
                        if (p.type != fTypeB) p.index = peptide_length(words[colInd[libPr]]) - p.index;
                        if (fragment_loss_info) {
                            auto &wl = words[colInd[libFrLoss]];
                            if (wl == "noloss") p.loss = loss_none;
                            else if (wl == "H2O") p.loss = loss_H2O;
                            else if (wl == "NH3") p.loss = loss_NH3;
                            else if (wl == "CO") p.loss = loss_CO;
                            else p.loss = loss_other;
                        } else p.loss = loss_none;
                    }
                }
                if (colInd[libFrExc] >= 0) {
                    auto &we = words[colInd[libFrExc]];
                    if (we.size()) if (we[0] == 'T' || we[0] == '1' || we[0] == 't') p.type |= fExclude;
                }

                bool decoy_fragment = false;
                if (colInd[libIsDecoy] >= 0) {
                    auto &w = words[colInd[libIsDecoy]];
                    if (w.length()) if (w[0] == '1' || w[0] == 't' || w[0] == 'T') decoy_fragment = true;
                }
                if (decoy_fragment && FastaSearch) continue; // do not use supplied decoys for the guide library
                if (decoy_fragment) gen_decoys = false;
                int charge = std::stoi(words[colInd[libCharge]]);
                id = to_canonical(words[colInd[libPr]], charge);
                if (id != prev_id || !prev_id.length()) {
                    ins = map.insert(std::pair<std::string, Entry>(id, e)).first;
                    prev_id = id;
                    ins->second.name = ins->first;

                    auto eg_name = (colInd[libEG] >= 0 ? words[colInd[libEG]] : to_eg(words[colInd[libPr]]));
                    auto egp = eg.insert(std::pair<std::string, int>(eg_name, eg.size()));
                    ins->second.eg_id = egp.first->second;

                    auto prot_id = (colInd[libPID] >= 0 ? words[colInd[libPID]] : "");
                    ins->second.prot = prot.insert(PG(prot_id)).first;
                    if (colInd[libPN] >= 0) ins->second.prot->names = words[colInd[libPN]];
                    if (colInd[libGenes] >= 0) ins->second.prot->genes = words[colInd[libGenes]];
                    if (colInd[libPT] >= 0 && words[colInd[libPT]].size()) {
                        char pt = words[colInd[libPT]][0];
                        if (pt == 'T' || pt == 't' || pt == '1' || pt == 'Y' || pt == 'y') ins->second.proteotypic = true;
                        else ins->second.proteotypic = false;
                    }
                }

                auto pep = !decoy_fragment ? (&(ins->second.target)) : (&(ins->second.decoy));
                pep->mz = to_double(words[colInd[libPrMz]]);
                pep->iRT = to_double(words[colInd[libiRT]]);
                if (colInd[libQ] >= 0) pep->lib_qvalue = to_double(words[colInd[libQ]]);
                pep->charge = charge;
                pep->fragments.push_back(p);
            }
            csv.close();

            eg.clear();
            entries.resize(map.size());
            precursors.resize(map.size());
            elution_groups.resize(map.size());
            i = 0;
            for (auto it = map.begin(); it != map.end(); it++, i++) {
                precursors[i] = it->first;
                entries[i] = it->second;
                entries[i].target.index = entries[i].decoy.index = i;
                if (entries[i].target.iRT < iRT_min) iRT_min = entries[i].target.iRT;
                if (entries[i].target.iRT > iRT_max) iRT_max = entries[i].target.iRT;
                it->second.prot->precursors.push_back(i);
                elution_groups[i] = entries[i].eg_id;
            }
            protein_ids.insert(protein_ids.begin(), prot.begin(), prot.end());

            for (i = 0; i < protein_ids.size(); i++)
                for (auto it = protein_ids[i].precursors.begin(); it != protein_ids[i].precursors.end(); it++)
                    entries[*it].pid_index = i;

            extract_proteins();

            if (!gen_decoys) {
                qDebug() << "WARNING: the library contains decoy precursors, DIA-NN's built-in decoy generation algorithm will be turned off (not recommended); performance is likely to be suboptimal\n";
                int dwt = 0, twd = 0;
                for (auto &e : entries) {
                    if (!e.target.fragments.size() && e.decoy.fragments.size()) dwt++;
                    if (e.target.fragments.size() && !e.decoy.fragments.size()) twd++;
                }
                if (dwt) qDebug() << "WARNING: " << dwt << " decoy precursors don't have matching target precursors; each decoy supplied in the library must correspond to a target with the same modified sequence and charge\n";
                if (twd) qDebug() << "WARNING: " << twd << " target precursors don't have matching decoy precursors; is this intended?\n";
            }
        }

        int ps = proteins.size(), pis = protein_ids.size();
        if (ps) if (!proteins[0].id.size()) ps--;
        if (pis) if (!protein_ids[0].ids.size()) pis--;

        elution_group_index();

        if (Verbose >= 1) Time(), qDebug() << "Spectral library loaded: "
                                        << ps << " protein isoforms, "
                                        << pis << " protein groups and "
                                        << entries.size() << " precursors in " << co_elution_index.size() << " elution groups.\n";

        return true;
    }

    void export_prosit() {
        prosit_file = remove_extension(out_lib_file) + std::string(".prosit.csv");
        std::ofstream out(prosit_file);
        if (out.fail()) {
            return;
        }
        out << "modified_sequence,collision_energy,precursor_charge\n";
        std::set<std::pair<std::pair<std::string, int>, float> > prs;
        for (auto &e : entries) prs.insert(std::pair<std::pair<std::string, int>, float>(std::pair<std::string, int>(to_prosit(e.name), e.target.charge), e.target.mz));
        for (auto &s : prs) out << s.first.first << ",30," << s.first.second << '\n';
        out.close();

    }

    void load_protein_annotations(Fasta &fasta) {
        std::set<PG> pid;
        std::string s;
        std::pair<std::string, std::vector<int> > key;
        for (int i = 0; i < entries.size(); i++) {
            bool prot_from_lib = false;
            if (OverwriteLibraryPGs || (entries[i].entry_flags & fFromFasta)) {
                auto stripped = get_aas(entries[i].name);
                if (stripped.size() < 5) prot_from_lib = true;
                else {
                    key.first = stripped.substr(0, 5);
                    auto pos = std::lower_bound(fasta.dict.begin(), fasta.dict.end(), key);
                    if (pos != fasta.dict.end()) if (pos->first == key.first) { // checks needed when !(entries[i].entry_flags & fFromFasta)
                            std::string ids;
                            for (auto &seq : pos->second) {
                                auto &sequence = fasta.sequences[seq];
                                for (int shift = 0; shift <= sequence.size() - stripped.size(); shift++) {
                                    int found = sequence.find(stripped, shift);
                                    if (found != std::string::npos) {
                                        bool l_cut = false, r_cut = false;
                                        shift = found;
                                        if (found == 0) l_cut = true;
                                        else if (found == 1 && NMetExcision && (sequence[0] == 'M' || sequence[0] == 'X')) l_cut = true;
                                        else l_cut = do_cut(sequence[found - 1], sequence[found]);
                                        if (l_cut) {
                                            if (found >= sequence.size() - stripped.size()) r_cut = true;
                                            else r_cut = do_cut(sequence[found + stripped.size() - 1], sequence[found + stripped.size()]);
                                        }
                                        if (l_cut && r_cut) {
                                            if (!ids.size()) ids = fasta.ids[seq];
                                            else if (!name_included(ids, fasta.ids[seq])) ids += std::string(";") + fasta.ids[seq];
                                            break;
                                        }
                                    } else break;
                                }
                            }
                            if (ids.size()) {
                                auto ins = pid.insert(PG(ids));
                                ins.first->precursors.push_back(i);
                            } else prot_from_lib = true;
                        } else prot_from_lib = true;
                }
            }
            if (prot_from_lib || (!(entries[i].entry_flags & fFromFasta) && !OverwriteLibraryPGs)) {
                auto ins = pid.insert(PG(protein_ids[entries[i].pid_index].ids));
                ins.first->precursors.push_back(i);
            }
        }

        protein_ids.clear();
        protein_ids.insert(protein_ids.begin(), pid.begin(), pid.end());
        pid.clear();

        for (int i = 0; i < protein_ids.size(); i++)
            for (auto it = protein_ids[i].precursors.begin(); it != protein_ids[i].precursors.end(); it++)
                entries[*it].pid_index = i;

        extract_proteins();
        if (Verbose >= 1) Time(), qDebug() << entries.size() << " precursors generated\n";

        if (ExportProsit) export_prosit();
    }

    void load(Fasta &fasta) {
        if (Verbose >= 1) Time(), qDebug() << "Processing FASTA\n";

        int i = entries.size(), cnt = 0;
        if (!i && !InSilicoRTPrediction) iRT_min = -1.0, iRT_max = 1.0;
        double default_iRT = (iRT_min + iRT_max) * 0.5;
        Entry e;
        e.lib = this; e.entry_flags = fFromFasta;
        for (auto &pep : fasta.peptides) {
            auto seq = get_sequence(pep);
            auto fragments = generate_fragments(seq, 1, loss_none, &cnt, MinFrAAs);
            auto aas = get_aas(pep);
            double mass = sum(seq) + proton + OH;
            int min_charge = Max(MinPrCharge, Max(1, (int)(mass / MaxPrMz)));
            int max_charge = Min(MaxPrCharge, Max(min_charge, Min(MaxGenCharge, (int)(mass / MinPrMz))));
            double iRT = InSilicoRTPrediction ? predict_irt(pep) : default_iRT;
            if (iRT < iRT_min) iRT_min = iRT;
            if (iRT > iRT_max) iRT_max = iRT;
            for (int charge = min_charge; charge <= max_charge; charge++) {
                double mz = proton + mass / (double)charge;
                if (mz < MinPrMz || mz > MaxPrMz) continue;
                e.target.init(mz, iRT, charge, i);
                e.target.fragments.clear();
                cnt = 0;
                for (auto &fr : fragments)
                    if (fr.mz >= MinFrMz && fr.mz <= MaxFrMz)
                        if (fr.type & fTypeY) cnt++;
                if (cnt >= MinGenFrNum) {
                    e.name = to_canonical(pep, charge);
                    auto pos = std::lower_bound(precursors.begin(), precursors.end(), e.name);
                    if (pos == precursors.end() || *pos != e.name) entries.push_back(e), i++;
                }
            }
        }

        assemble_elution_groups();
        elution_group_index();

        load_protein_annotations(fasta);
    }

    void library_protein_stats() {
        int ns = names.size(), gs = genes.size();
        if (ns) if (!names[0].size()) {
                ns--;
                if (Verbose >= 1) Time(), qDebug() << "Protein names missing for some isoforms\n";
            }
        if (gs) {
            if (!genes[0].size()) gs--;
            if (Verbose >= 1) Time(), qDebug() << "Gene names missing for some isoforms\n";
        }
        if (Verbose >= 1) Time(), qDebug() << "Library contains " << ns << " proteins, and " << gs << " genes\n";
    }

    void annotate() {
        std::set<std::string> name(names.begin(), names.end()), gene(genes.begin(), genes.end());
        for (auto &p : proteins) name.insert(p.name), gene.insert(p.gene);
        names.clear(), genes.clear();
        names.insert(names.begin(), name.begin(), name.end());
        genes.insert(genes.begin(), gene.begin(), gene.end());
        for (auto &p : proteins) {
            p.name_index = std::distance(names.begin(), std::lower_bound(names.begin(), names.end(), p.name));
            p.gene_index = std::distance(genes.begin(), std::lower_bound(genes.begin(), genes.end(), p.gene));
        }
        library_protein_stats();
    }

    void annotate_pgs(std::vector<PG> &pgs) {
        for (auto &pid : pgs) pid.annotate(proteins, names, genes);
    }

    void save(const std::string &file_name, std::vector<int> * ref, bool searched, bool decoys = false) {

        std::ofstream out(file_name, std::ofstream::out);
        if (out.fail()) { qDebug() << "ERROR: cannot write to " << QString::fromStdString(file_name) << ". Check if the destination folder is write-protected or the file is in use\n"; return; }
        out << "FileName\tPrecursorMz\tProductMz\tTr_recalibrated\ttransition_name\tLibraryIntensity\ttransition_group_id\tdecoy\tPeptideSequence\tProteotypic\tQValue\tPGQValue\t";
        out << "ProteinGroup\tProteinName\tGenes\tFullUniModPeptideName\tModifiedPeptide\t";
        out << "PrecursorCharge\tPeptideGroupLabel\tUniprotID\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\tFragmentLossType\tExcludeFromAssay\n";
        out.precision(8);

        auto &prot = (InferPGs && searched) ? protein_groups : protein_ids;
        if (!protein_ids.size()) protein_ids.resize(1);
        if (!prot.size()) prot.resize(1);

        bool rec_info = false;
        int cnt = -1, skipped = 0, empty = 0, targets_written = 0, decoys_written = 0;
        for (auto &it : entries) {
            cnt++;
            if (searched && it.best_run < 0 && ((FastaSearch && (it.entry_flags & fFromFasta)) || !SaveOriginalLib)) continue;
            if (!searched && it.target.fragments.size() < MinF) continue;
            if (ref != NULL) {
                auto pos = std::lower_bound(ref->begin(), ref->end(), cnt);
                if (pos == ref->end()) continue;
                if (*pos != cnt) continue;
            }
            auto seq = get_aas(it.name);
            float gen_acc = 0;

            bool recognise = false;
            std::vector<Ion> ions;
            for (auto &f : it.target.fragments) if (!(f.type & 3) || f.charge <= 0) {
                    recognise = true;
                    if (!rec_info) {
                        rec_info = true;
                        if (Verbose >= 1) Time(), qDebug() << "Some precursors lack fragment annotation; annotating...\n";
                    }
                    break;
                }
            if (recognise) {
                ions = recognise_fragments(it.name, it.target.fragments, gen_acc, it.target.charge, ExportRecFragOnly, MaxRecCharge, MaxRecLoss);
                if (ExportRecFragOnly && !searched) {
                    int rec = 0;
                    for (int fr = 0; fr < ions.size(); fr++) {
                        ions[fr].type = it.target.fragments[fr].type;
                        if (ions[fr].charge) rec++;
                    }
                    if (rec < MinF) continue;
                }
            } else {
                auto &frs = it.target.fragments;
                ions.resize(frs.size());
                for (int i = 0; i < ions.size(); i++) ions[i].init(frs[i]), ions[i].type = frs[i].type;
            }
            auto &pep = it.target;
            if (!pep.fragments.size()) {
                empty++;
                continue;
            }
            auto pep_name = to_canonical(it.name);
            auto name = pep_name + std::to_string(pep.charge);
            for (int dc = 0; dc <= 1; dc++) {
                if (!decoys && dc) continue;
                std::string prefix = dc ? std::string("DECOY") : "";
                auto &pep = dc ? it.decoy : it.target;
                int fr_cnt = 0;
                for (int fr = 0; fr < pep.fragments.size(); fr++) {
                    int fr_type, fr_loss, fr_num, fr_charge = ions[fr].charge;
                    float fr_mz = pep.fragments[fr].mz;
                    if (fr_charge) {
                        fr_type = ions[fr].type;
                        fr_loss = ions[fr].loss;
                        fr_num = ions[fr].index;
                        if (ExportRecFragOnly && !dc) fr_mz = ions[fr].mz;
                    } else { // unrecognised
                        fr_type = fr_loss = fr_num = 0;
                        if (ExportRecFragOnly) continue;
                    }
                    if (!(it.entry_flags & fFromFasta)) {
                        bool skip = false;
                        for (auto pos = ions.begin(); pos < ions.begin() + fr; pos++)
                            if (pos->index == fr_num && pos->type == fr_type && pos->charge == fr_charge && pos->loss == fr_loss) {
                                skip = true;
                                break;
                            }
                        if (skip) continue;
                    }
                    if (pep.fragments[fr].height < MinGenFrInt) continue;
                    int pg = (InferPGs && searched) ? it.pg_index : it.pid_index;
                    auto fr_name = name + std::string("_") + std::to_string(char_from_type[fr_type & 3])
                                   + std::string("_") + std::to_string(fr_charge) + std::string("_")
                                   + std::to_string(fr_loss) + std::string("_") + std::to_string(fr_num);
                    out << ((it.best_run >= 0) ? ms_files[it.best_run].c_str() : lib_file.c_str()) << '\t'
                        << pep.mz << '\t'
                        << fr_mz << '\t'
                        << pep.iRT << '\t'
                        << (prefix + fr_name).c_str() << '\t'
                        << pep.fragments[fr].height << '\t'
                        << (prefix ).c_str() << '\t'
                        << dc << '\t'
                        << seq << '\t'
                        << (int)it.proteotypic << '\t'
                        << it.qvalue << '\t'
                        << it.pg_qvalue << '\t'
                        << (prefix + prot[pg].ids).c_str() << '\t'
                        << prot[pg].names.c_str() << '\t'
                        << prot[pg].genes.c_str() << '\t'
                        << pep_name.c_str() << '\t'
                        << pep_name.c_str() << '\t'
                        << pep.charge << '\t'
                        << (prefix + pep_name).c_str() << '\t'
                        << protein_ids[it.pid_index].ids.c_str() << '\t'
                        << char_from_type[fr_type & 3] << '\t'
                        << fr_charge << '\t'
                        << ((fr_type & fTypeB) ? fr_num : seq.size() - fr_num) << '\t'
                        << name_from_loss[fr_loss].c_str() << "\t"
                        << ((fr_type & fExclude) ? "True" : "False") << "\n";
                    fr_cnt++;
                }
                if (!fr_cnt) empty++;
                else {
                    if (!dc) targets_written++;
                    else decoys_written++;
                }
            }
        }

        out.close();
        if (Verbose >= 1) {
            if (!decoys_written) Time(), qDebug() << targets_written << " precursors saved\n";
            else Time(), qDebug() << targets_written << " target and " << decoys_written << " decoy precursors saved\n";
        }
        if (skipped) qDebug() << "WARNING: " << skipped << " precursors with unrecognised library fragments were skipped\n";
        if (empty) qDebug() << "WARNING: " << empty << " precursors without any fragments annotated were skipped\n";
    }

    bool save(const std::string &file_name) {
        if (Verbose >= 1) Time(), qDebug() << "Saving the library to " << QString::fromStdString(file_name) << "\n";
        std::ofstream speclib(file_name, std::ofstream::binary);
        if (speclib.fail()) {
            qDebug() << "Could not save " << QString::fromStdString(file_name) << "\n";
            return false;
        } else {
            write(speclib);
            speclib.close();
            return true;
        }
    }

    void init(bool decoys) {
        bool corrected_charge = false;
        for (auto it = entries.begin(); it != entries.end(); it++) if (it->lock.set()) {
                if (it->target.charge <= 0) it->target.charge = 2, corrected_charge = true;
                if (gen_decoys && decoys) it->generate_decoy();
                if (it->decoy.charge <= 0 && decoys) it->decoy.charge = 2, corrected_charge = true;
                else if (UseIsotopes && gen_charges && !(it->entry_flags & fFromFasta)) it->annotate_charges();
                it->init();
            }
        if (corrected_charge) qDebug() << "WARNING: at least one library precursor had non-positive charge; corrected to charge = 2\n";
    }

    void initialise(bool decoys) {
        int i;

        if (Predictor) {
#ifdef PREDICTOR
            generate_predictions();
			predictor_file = remove_extension(out_lib_file) + std::string(".predicted.speclib");

			// save the full list of precursors to the .predicted.speclib file
			auto old_precursors = precursors;
			precursors.resize(entries.size());
			for (i = 0; i < entries.size(); i++) precursors[i] = entries[i].name;

			PredictorSaved = save(predictor_file);
			precursors = old_precursors;
#endif
        }
        skipped = 0;
        if (Verbose >= 1) Time(), qDebug() << "Initialising library\n";

        init(decoys);

        for (i = 0; i < entries.size(); i++) entries[i].lock.free();
//        if (skipped) qDebug() << "WARNING: " << skipped
//                           << " library precursors were skipped due to unrecognised fragments when generating decoys.\n";
//        if (ExportLibrary) save(out_lib_file, NULL, false, ExportDecoys);
    }

    void infer_proteins() { // filters precursors and proteins at q_cutoff level and infers protein groups
        int i;
        if (PGsInferred || !protein_ids.size() || !proteins.size()) return;
        protein_groups.clear();

        if (Verbose >= 1) Time(), qDebug() << "Assembling protein groups\n";

        // get IDs
        std::vector<std::pair<int, int> > IDs;
        std::vector<std::vector<int> > sets(proteins.size());
        std::vector<int> id_class;
        std::vector<int> max_protein_score(proteins.size(), (int)-INF);

        for (auto it = info.map.begin(); it != info.map.end(); it++) {
            auto v = &(it->second);

            for (auto jt = v->begin(); jt != v->end(); jt++) {
                int s = jt->first;
                auto pr = &(jt->second.pr);
                if (pr->qvalue <= ProteinIDQvalue) {
                    auto &pid = protein_ids[entries[pr->index].pid_index];

                    bool first = true;
                    id_class.clear(); id_class.resize(pid.proteins.size(), 0);
                    int i = -1, max_class = -1;
                    for (auto &p : pid.proteins) {
                        i++;
                        if (SwissProtPriority) if (proteins[p].swissprot) id_class[i] += 1;
                        auto pos = std::lower_bound(info.proteins[s].begin(), info.proteins[s].end(), p);
                        if (pos != info.proteins[s].end()) if (*pos == p) id_class[i] += 10;
                        if (PGLevel == 0) if (!proteins[p].id.size()) id_class[i] -= 100;
                        if (PGLevel == 1) if (!proteins[p].name.size()) id_class[i] -= 100;
                        if (PGLevel == 2) if (!proteins[p].gene.size()) id_class[i] -= 100;
                        if (id_class[i] > max_class) max_class = id_class[i];
                        if (id_class[i] > max_protein_score[p]) max_protein_score[p] = id_class[i];
                    }

                    i = -1;
                    if (max_class >= 0) for (auto &p : pid.proteins) {
                            i++;
                            if (id_class[i] < max_class) continue;

                            if (first) {
                                sets[p].push_back(IDs.size());
                                IDs.push_back(std::pair<int, int>(pr->index, s));
                            } else sets[p].push_back(IDs.size() - 1);
                            first = false;
                        }
                }
            }
        }

        auto groups = bipartite_set_cover(sets, IDs.size());

        // for each precursor merge all groups associated with it
        std::vector<std::set<int> > pgs(entries.size());
        for (auto &g : groups) for (auto &e : g.e) {
                auto &pr = pgs[IDs[e].first];
                for (auto &s : g.s)
                    pr.insert(s);
            }

        // remove proteins for which there is less evidence than for other proteins in the same group
        for (auto &pr : pgs) {
            int max = -INF;
            for (auto &prot : pr) if (max_protein_score[prot] > max) max = max_protein_score[prot];
            for (auto &prot : pr) if (max_protein_score[prot] < max) pr.erase(prot);
        }

        std::set<PG> pg;
        std::string name;
        for (i = 0; i < entries.size(); i++) {
            auto &pid = protein_ids[entries[i].pid_index];
            if (!pgs[i].size()) {
                if (!pid.proteins.size()) name = pid.ids;
                else {
                    auto it = pid.proteins.begin();
                    name = proteins[*it].id;
                    for (it++; it != pid.proteins.end(); it++)
                        name += ';' + proteins[*it].id;
                }
            } else {
                auto it = pgs[i].begin();
                name = proteins[*it].id;
                for (it++; it != pgs[i].end(); it++)
                    name += ';' + proteins[*it].id;
            }
            auto ins = pg.insert(PG(name));
            if (!pgs[i].size() && ins.second) ins.first->proteins = pid.proteins;
            for (auto &p : pgs[i]) ins.first->proteins.insert(p);
            ins.first->precursors.push_back(i);
        }
        protein_groups.clear();
        protein_groups.insert(protein_groups.begin(), pg.begin(), pg.end());

        for (i = 0; i < protein_groups.size(); i++)
            for (auto it = protein_groups[i].precursors.begin(); it != protein_groups[i].precursors.end(); it++)
                entries[*it].pg_index = i;

        bool annotate = false;
        if (names.size()) if (names[names.size() - 1] != "" && names[names.size() - 1] != " ") annotate = true;
        if (genes.size()) if (genes[genes.size() - 1] != "" && genes[genes.size() - 1] != " ") annotate = true;
        if (annotate) annotate_pgs(protein_groups);
    }

    int generate_gene_groups() {
        std::set<std::string> gg;
        std::vector<std::string> ggs;
        gg_index.clear(); gg_index.resize(entries.size(), 0);
        auto &prot = (InferPGs ? protein_groups : protein_ids);

        for (auto &pg : prot) gg.insert(pg.genes);
        ggs.clear();  ggs.insert(ggs.begin(), gg.begin(), gg.end()); gg.clear();
        if (!ggs.size()) ggs.resize(1);
        int n = ggs.size();

        for (auto &pg : prot) {
            int ind = std::distance(ggs.begin(), std::lower_bound(ggs.begin(), ggs.end(), pg.genes));
            for (auto &pr : pg.precursors) gg_index[pr] = ind;
        }

        gene_groups.clear(), gene_groups.resize(n);
        for (int i = 0; i < n; i++) gene_groups[i].ids = ggs[i];
        for (int i = 0; i < gg_index.size(); i++) gene_groups[gg_index[i]].precursors.push_back(i);
        return n;
    }

    void group_qvalues(int stage = 0) {

        auto &prot = stage ? gene_groups : (InferPGs ? protein_groups : protein_ids);
        int i, s, n_s = info.n_s, n_p = prot.size();
        if (!n_p && Verbose >= 1) {
            Time(), qDebug() << "No protein annotation, skipping protein q-value calculation\n";
            for (auto it = info.map.begin(); it != info.map.end(); it++) {
                auto v = &(it->second);
                for (auto jt = v->begin(); jt != v->end(); jt++) {
                    auto pr = &(jt->second.pr);
                    if (stage == 0) pr->pg_qvalue = 1.0;
                    else pr->gg_qvalue = 1.0;
                }
            }
            if (stage < 1) group_qvalues(stage + 1);
            return;
        }
        std::vector<float> ts(n_s * n_p, 0.0), ds(n_s * n_p, 0.0), tq(n_s * n_p, 0.0), dq(n_s * n_p, 0.0),
                qvalue[2], targets(n_p), decoys(n_p);
        std::vector<int> prs(n_s * n_p, 0);
        if (Verbose >= 1 && stage == 0) Time(), qDebug() << "Calculating q-values for protein and gene groups\n";

        for (i = 0; i < entries.size(); i++) {
            int p = stage ? gg_index[i] : (InferPGs ? entries[i].pg_index : entries[i].pid_index);
            int ind = i / sizeof(int);
            int shift = i % sizeof(int);
            for (s = 0; s < n_s; s++) if (info.detectable_precursors[s * info.n_det_pr_block + ind] & (1 << shift)) prs[s * n_p + p]++;
        }

        for (auto it = info.map.begin(); it != info.map.end(); it++) {
            auto v = &(it->second);

            for (auto jt = v->begin(); jt != v->end(); jt++) {
                int s = jt->first;
                auto pr = &(jt->second.pr);

                int p = stage ? gg_index[pr->index] : (InferPGs ? entries[pr->index].pg_index : entries[pr->index].pid_index);
                if (!prot[p].ids.size()) continue;

                if (prs[s * n_p + p] < E) qDebug() << "ERROR: undetectable precursor id reported in the .quant file\n";
                float err = pr->lib_qvalue + Max(pr->lib_qvalue, pr->dratio);
                float sc = Max(0.0, 1.0 - err);
                if (pr->qvalue <= ProteinQCalcQvalue) {
                    if (sc > tq[s * n_p + p]) tq[s * n_p + p] = sc;
                    ts[s * n_p + p] -= log(Max(E, Min(1.0, err * prs[s * n_p + p])));
                }
            }
        }

        for (auto it = info.decoy_map.begin(); it != info.decoy_map.end(); it++) {
            auto v = &(it->second);

            for (auto jt = v->begin(); jt != v->end(); jt++) {
                int s = jt->first;
                auto pr = &(jt->second);

                int p = stage ? gg_index[pr->index] : (InferPGs ? entries[pr->index].pg_index : entries[pr->index].pid_index);
                if (!prot[p].ids.size()) continue;

                if (prs[s * n_p + p] < E) qDebug() << "ERROR: undetectable precursor id reported in the .quant file\n";
                float err = pr->lib_qvalue + Max(pr->lib_qvalue, pr->dratio);
                float sc = Max(0.0, 1.0 - err);
                if (pr->qvalue <= ProteinQCalcQvalue) {
                    if (sc > dq[s * n_p + p]) dq[s * n_p + p] = sc;
                    ds[s * n_p + p] -= log(Max(E, Min(1.0, err * prs[s * n_p + p])));
                }
            }
        }

        qvalue[0].resize(n_s * n_p, 1.0);
        qvalue[1].resize(n_s * n_p, 1.0);
        for (s = 0; s < n_s; s++) {
            for (int at = 0; at <= 1; at++) {
                auto tsc = (at == 0 ? &ts[s * n_p] : &tq[s * n_p]);
                auto dsc = (at == 0 ? &ds[s * n_p] : &dq[s * n_p]);
                auto q = &qvalue[at][s * n_p];
                memcpy(&targets[0], tsc, n_p * sizeof(float));
                memcpy(&decoys[0], dsc, n_p * sizeof(float));

                std::sort(targets.begin(), targets.end());
                std::sort(decoys.begin(), decoys.end());

                std::map<float, float> cs_qv;
                for (i = 0; i < n_p; i++) {
                    if (tsc[i] < E) continue;
                    float sc = tsc[i];
                    auto pD = std::lower_bound(decoys.begin(), decoys.end(), sc);
                    if (pD > decoys.begin()) if (*(pD - 1) > E) sc = *(--pD);
                    auto pT = std::lower_bound(targets.begin(), targets.end(), sc);

                    q[i] = Min(1.0, ((double)std::distance(pD, decoys.end())) / Max(1.0, (double)std::distance(pT, targets.end())));

                    auto pair = std::pair<float, float>(tsc[i], q[i]);
                    auto pos = cs_qv.insert(pair);
                    if (pos.second) {
                        if (pos.first != cs_qv.begin() && std::prev(pos.first)->second < pair.second) pos.first->second = std::prev(pos.first)->second;
                        else for (auto jt = std::next(pos.first); jt != cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
                    }
                }

                for (i = 0; i < n_p; i++) {
                    if (q[i] < 1.0) {
                        auto pos = cs_qv.lower_bound(tsc[i]);
                        q[i] = pos->second;
                    }
                }
            }
            auto q0 = &qvalue[0][s * n_p], q1 = &qvalue[1][s * n_p];
            for (i = 0; i < n_p; i++) if (q1[i] < q0[i]) q0[i] = q1[i];
        }

        for (auto it = info.map.begin(); it != info.map.end(); it++) {
            auto v = &(it->second);

            for (auto jt = v->begin(); jt != v->end(); jt++) {
                int s = jt->first;
                auto pr = &(jt->second.pr);
                int p = stage ? gg_index[pr->index] : (InferPGs ? entries[pr->index].pg_index : entries[pr->index].pid_index);
                float qv = qvalue[0][s * n_p + p];
                if (stage == 0) pr->pg_qvalue = qv;
                else pr->gg_qvalue = qv;
            }
        }

        if (stage < 1) group_qvalues(stage + 1);
    }

    void quantify_proteins(int N, double q_cutoff, int stage = 0) { // top N method for protein quantification
        if (InferPGs && !stage) infer_proteins();
        auto &prot = (InferPGs ? protein_groups : protein_ids);
        int i, j, k, pass, n = stage ? genes.size() : prot.size();

        if (Verbose >= 1 && stage == 0) Time(), qDebug() << "Quantifying proteins\n";

        if (stage == 1) n = generate_gene_groups();
        if (!n) {
            if (stage < 1) quantify_proteins(N, q_cutoff, stage + 1);
            return;
        }
        auto &proteins = stage ? gene_groups : prot;

        std::vector<float> max_q(n * info.n_s, 1.0), top_l(n * info.n_s * N, 0.0), quant(n * info.n_s, 0.0), norm(n * info.n_s, 0.0);

        for (pass = 0; pass <= 3; pass++) {
            for (auto it = info.map.begin(); it != info.map.end(); it++) {
                auto v = &(it->second);

                for (auto jt = v->begin(); jt != v->end(); jt++) {
                    int s = jt->first;
                    auto pr = &(jt->second.pr);

                    if (pass == 0) pr->max_lfq = pr->max_lfq_unique = 0.0;
                    int ggi = (stage == 1) ? gg_index[pr->index] : 0;
                    if (stage == 1) if (!gene_groups[gg_index[pr->index]].ids.size()) {
                            if (pass == 0) pr->gene_quantity = pr->gene_norm = 0.0;
                            continue;
                        }
                    int pg = InferPGs ? entries[pr->index].pg_index : entries[pr->index].pid_index;
                    int gi = 0;

                    float quantity = pr->norm;
                    float q = pr->qvalue;

                    int index = (stage == 1 ? ggi : pg) * info.n_s + s;
                    auto pos_q = &(max_q[index]);
                    auto pos_l = &(top_l[index * N]);

                    if (pass == 0) {
                        if (q < *pos_q && *pos_q > q_cutoff) *pos_q = Max(q_cutoff, q);
                    } else if (pass == 1) {
                        if (q <= *pos_q) for (i = 0; i < N; i++) if (quantity > pos_l[i]) {
                                    for (j = N - 1; j > i; j--) pos_l[j] = pos_l[j - 1];
                                    pos_l[i] = quantity;
                                    break;
                                }
                    } else if (pass == 2) {
                        if (q <= *pos_q && quantity >= pos_l[N - 1]) quant[index] += pr->quantity, norm[index] += pr->norm;
                    } else if (pass == 3) {
                        if (stage == 1) pr->gene_quantity = quant[index], pr->gene_norm = norm[index];
                        else pr->pg_quantity = quant[index], pr->pg_norm = norm[index];
                    }
                }
            }
        }

        if (stage == 1 && MaxLFQ) { // gene groups
            if (Verbose >= 1 && ms_files.size() > 2000)
                Time(), qDebug() << "MaxLFQ protein quantification can be slow for very large experiments; if it's taking too long, consider stopping DIA-NN and rerunning with --no-maxlfq and --use-quant options - the latter to reuse the existing .quant files\n";

            std::vector<Lock> locks(proteins.size());
            std::vector<std::thread> thr;

        }

        if (stage == 1 && MaxLFQ) { // unique genes
            std::vector<Lock> locks(proteins.size());
            std::vector<std::thread> thr;

        }

        if (stage < 1) quantify_proteins(N, q_cutoff, stage + 1);
        else group_qvalues();
    }

    void report_pr_matrix(const std::string &file_name) { // must only be called after the regular report() below
        if (Verbose >= 1) Time(), qDebug() << "Saving precursor levels matrix\n";
        std::ofstream out(file_name, std::ofstream::out);

        out << "Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tProteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id";
        for (auto &f : ms_files) out << '\t' << f; out << '\n';

        auto &prot = InferPGs ? protein_groups : protein_ids;
        std::vector<std::pair<int, float> > quants;

        for (auto it = info.map.begin(); it != info.map.end(); it++) {
            auto v = &(it->second);
            auto entry = &(entries[it->first]);
            int pg = InferPGs ? entry->pg_index : entry->pid_index;

            quants.clear();
            for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                double q = Max(jt->second.pr.qvalue, jt->second.pr.pg_qvalue);
                if (q <= MatrixQValue) quants.push_back(std::pair<int, float>(jt->first, jt->second.pr.norm));
            }
            if (!quants.size()) continue;

            out << prot[pg].ids.c_str() << '\t'
                << protein_ids[entry->pid_index].ids.c_str() << '\t'
                << prot[pg].names.c_str() << '\t'
                << prot[pg].genes.c_str() << '\t';
            if (prot[pg].proteins.size()) out << proteins[*(prot[pg].proteins.begin())].description.c_str() << '\t';
            else out << '\t';
            out << (int)entry->proteotypic << '\t'
                << get_aas(entry->name).c_str() << '\t'
                << pep_name(entry->name).c_str() << '\t'
                << entry->target.charge << '\t'
                << entry->name.c_str();

            int pos = 0;
            std::sort(quants.begin(), quants.end());
            for (auto &q : quants) {
                for (; pos < q.first; pos++) out << '\t'; pos++;
                out << '\t' << q.second;
            }
            for (; pos < ms_files.size(); pos++) out << '\t';
            out << '\n';
        }

        out.close();


    }

    void report_pg_matrix(const std::string &file_name, bool genes = false, bool unique = false) { // must only be called after the regular report() below
        if (Verbose >= 1) {
            if (!genes) Time(), qDebug() << "Saving protein group levels matrix\n";
            else if (!unique) Time(), qDebug() << "Saving gene group levels matrix\n";
            else Time(), qDebug() << "Saving unique genes levels matrix\n";
        }
        std::ofstream out(file_name, std::ofstream::out);

        out << "Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description";
        for (auto &f : ms_files) out << '\t' << f; out << '\n';
        if (ms_files.size() != info.n_s) {
            qDebug() << "ERROR: algorithmic error: info.n_s == " << info.n_s << ", while ms_files.size() == " << ms_files.size() << "\n";
            assert(ms_files.size() == info.n_s);
        }

        double q, qt;
        auto &prot = InferPGs ? protein_groups : protein_ids;
        std::vector<std::pair<int, float> > quants(ms_files.size());

        int i, j, n_pg = 0;
        for (auto it = info.map.begin(); it != info.map.end(); it++) {
            auto entry = &(entries[it->first]);
            int pg = InferPGs ? entry->pg_index : entry->pid_index;
            if (pg >= n_pg) n_pg = pg + 1;
        }
        if (!n_pg) {
            out.close();
            return;
        }
        std::vector<float> PG(n_pg * info.n_s, 0.0);
        std::vector<int> pid_index(n_pg, 0);

        for (auto it = info.map.begin(); it != info.map.end(); it++) {
            auto v = &(it->second);
            auto entry = &(entries[it->first]);
            int pg = InferPGs ? entry->pg_index : entry->pid_index;
            pid_index[pg] = entry->pid_index;
            int shift = pg * info.n_s;

            if (genes && unique) if (!entry->proteotypic) continue;

            for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
                if (PG[shift + jt->first] > E) continue;
                if (!genes) q = Max(jt->second.pr.qvalue, jt->second.pr.pg_qvalue), qt = jt->second.pr.pg_norm;
                else if (!unique) q = Max(jt->second.pr.qvalue, jt->second.pr.gg_qvalue), qt = jt->second.pr.max_lfq;
                else q = Max(jt->second.pr.qvalue, jt->second.pr.protein_qvalue), qt = jt->second.pr.max_lfq_unique;
                if (q <= MatrixQValue) PG[shift + jt->first] = qt;
            }
        }

        for (i = 0; i < n_pg; i++) {
            bool found = false;
            int shift = i * info.n_s;
            for (j = 0; j < info.n_s; j++) if (PG[shift + j] > E) {
                    found = true;
                    break;
                }
            if (!found) continue;

            out << prot[i].ids.c_str() << '\t'
                << protein_ids[pid_index[i]].ids.c_str() << '\t'
                << prot[i].names.c_str() << '\t'
                << prot[i].genes.c_str();
            if (prot[i].proteins.size()) out << '\t' << proteins[*(prot[i].proteins.begin())].description.c_str();
            else out << '\t';

            for (j = 0; j < info.n_s; j++) {
                if (PG[shift + j] > E) out << '\t' << PG[shift + j];
                else out << '\t';
            }
            out << '\n';
        }

        out.close();

    }


#if (HASH > 0)
    unsigned int hash() {
		unsigned int res = 0;
		for (int i = 0; i < entries.size(); i++) res ^= entries[i].hash();
		return res;
	}
#endif
};

class Profile {
public:
    std::vector<QuantEntry> entries;
    std::vector<float> tq, dq;

    void gen_decoy_list() {
        std::sort(dq.begin(), dq.end());
        auto pos = std::lower_bound(dq.begin(), dq.end(), ReportQValue);
        dq.resize(std::distance(dq.begin(), pos));
        std::sort(tq.begin(), tq.end());
        pos = std::lower_bound(tq.begin(), tq.end(), ReportQValue);
        tq.resize(std::distance(tq.begin(), pos));

        std::map<float, float> cs_qv;

        for (int i = 0; i < tq.size(); i++) {
            int n_targets = i + 1;
            int n_decoys = std::distance(dq.begin(), std::lower_bound(dq.begin(), dq.end(), tq[i]));
            float q = Min(1.0, ((double)Max(1, n_decoys)) / (double)Max(1, n_targets));

            auto pair = std::pair<float, float>(-tq[i], q);
            auto pos = cs_qv.insert(pair);
            if (pos.second) {
                if (pos.first != cs_qv.begin() && std::prev(pos.first)->second < pair.second) pos.first->second = std::prev(pos.first)->second;
                else for (auto jt = std::next(pos.first); jt != cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
            }
        }

        for (auto &e : entries) {
            if (e.index < 0) continue;
            auto pos = cs_qv.lower_bound(-(e.pr.qvalue + e.pr.lib_qvalue));
            if (pos == cs_qv.end()) e.pr.profile_qvalue = 1.0;
            else e.pr.profile_qvalue = pos->second;
        }
    }

    Profile(std::vector<std::string> &files, int lib_size = 0) {
        entries.resize(MaxLibSize);

        bool decoy_list = FastaSearch || GenSpecLib;
        if (decoy_list) dq.resize(MaxLibSize, 1.0), tq.resize(MaxLibSize, 1.0);

        for (int i = 0; i < files.size(); i++) {
            Quant Qt;
            auto Q = &Qt;
            if (!QuantInMem) {
                std::ifstream in(files[i] + std::string(".quant"), std::ifstream::binary);
                if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(files[i]) + std::string(".quant"), std::ifstream::binary);
                Q->read(in, lib_size);
                in.close();
            } else Q = &(quants[i]);

            for (auto it = Q->entries.begin(); it != Q->entries.end(); it++) {
                auto pos = it->pr.index;
                if (pos >= entries.size()) entries.resize(pos + 1 + pos / 2);
                if (entries[pos].index < 0) entries[pos] = *it;
                else if (it->pr.qvalue < entries[pos].pr.qvalue) entries[pos].pr = it->pr;

                if (decoy_list) {
                    float sc = Max(it->pr.qvalue, it->pr.lib_qvalue) + it->pr.lib_qvalue;
                    if (pos >= tq.size()) tq.resize(pos + 1 + pos / 2, 1.0);
                    if (sc <= ReportQValue && sc < tq[pos]) tq[pos] = sc;
                }
            }
            if (decoy_list) {
                for (auto it = Q->decoys.begin(); it != Q->decoys.end(); it++) {
                    auto pos = it->index;
                    float sc = Max(it->qvalue, it->lib_qvalue) + it->lib_qvalue;
                    if (pos >= dq.size()) dq.resize(pos + 1 + pos / 2, 1.0);
                    if (sc <= ReportQValue && sc < dq[pos]) dq[pos] = sc;
                }
            }
            if (!QuantInMem) std::vector<QuantEntry>().swap(Q->entries), std::vector<DecoyEntry>().swap(Q->decoys);
        }

        if (decoy_list) gen_decoy_list();
        std::vector<float>().swap(tq); std::vector<float>().swap(dq);
    }

    Profile(Library * lib) {
        entries.resize(MaxLibSize);

        bool decoy_list = FastaSearch || GenSpecLib;
        if (decoy_list) dq.resize(MaxLibSize, 1.0), tq.resize(MaxLibSize, 1.0);

        auto &info = lib->info;
        for (auto it = info.map.begin(); it != info.map.end(); it++) {
            auto v = &(it->second);

            for (auto jt = v->begin(); jt != v->end(); jt++) {
                int s = jt->first;
                auto &qr = jt->second;
                auto pos = qr.pr.index;
                if (pos >= entries.size()) entries.resize(pos + 1 + pos / 2);
                if (entries[pos].index < 0) entries[pos] = qr;
                else if (qr.pr.qvalue < entries[pos].pr.qvalue) entries[pos] = qr;

                if (decoy_list) {
                    float sc = Max(qr.pr.qvalue, qr.pr.lib_qvalue) + qr.pr.lib_qvalue;
                    if (pos >= tq.size()) tq.resize(pos + 1 + pos / 2, 1.0);
                    if (sc <= ReportQValue && sc < tq[pos]) tq[pos] = sc;
                }
            }
        }

        if (decoy_list) {
            for (auto it = info.decoy_map.begin(); it != info.decoy_map.end(); it++) {
                auto v = &(it->second);

                for (auto jt = v->begin(); jt != v->end(); jt++) {
                    int s = jt->first;
                    auto pr = &(jt->second);
                    int pos = pr->index;
                    float sc = Max(pr->qvalue, pr->lib_qvalue) + pr->lib_qvalue;
                    if (pos >= dq.size()) dq.resize(pos + 1 + pos / 2, 1.0);
                    if (sc <= ReportQValue && sc < dq[pos]) dq[pos] = sc;
                }
            }
        }

        if (decoy_list) gen_decoy_list();
        std::vector<float>().swap(tq); std::vector<float>().swap(dq);
    }
};

void annotate_library(Library &lib, Fasta &fasta) {
    for (auto &e : lib.proteins) {
        auto pos = std::lower_bound(fasta.proteins.begin(), fasta.proteins.end(), e);
        if (pos != fasta.proteins.end()) if (pos->id == e.id) e.name = pos->name, e.gene = pos->gene, e.description = pos->description, e.swissprot = pos->swissprot;
    }
    lib.annotate();
    lib.annotate_pgs(lib.protein_ids);
    for (auto &le : lib.entries) {
        if (lib.protein_ids[le.pid_index].proteins.size() == 1) le.proteotypic = true;
        else {
            if (PGLevel == 0) le.proteotypic = false;
            else if (PGLevel == 1) le.proteotypic = (lib.protein_ids[le.pid_index].name_indices.size() == 1);
            else if (PGLevel == 2) le.proteotypic = (lib.protein_ids[le.pid_index].gene_indices.size() == 1);
        }
    }
}

const int fFound = 1;
const int fDecoy = 1 << 1;
const int fInit = 1 << 2;
const int fTranslated = 1 << 3;


#endif //PYTHIADIACPP_MSFRAME_H