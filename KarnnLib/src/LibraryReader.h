//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_LIBRARYREADER_H
#define PYTHIADIACPP_LIBRARYREADER_H

#include "KarnnLib_Exports.h"

#include "ReadWriteCommonFunctons.h"

#include <map>
#include <set>
#include <string>
#include <vector>

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

#include "cpu_info.h"
#include "Ion.h"
#include "Isoform.h"
#include "LibraryCommon.h"
#include "Modifications.h"
#include "Peptide.h"
#include "Product.h"
#include "ProteinGroup.h"


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

class KARNNLIB_EXPORTS LibraryReader {

};

class KARNNLIB_EXPORTS Library {

    using PG = ProteinGroup;

public:

    std::string name;
    std::string fasta_names;
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

            const auto peptideSortLogicHeightAsc = [](const Product &l, const Product &r) {return l.height < r.height;};
            std::sort(target.fragments.rbegin(), target.fragments.rend(), peptideSortLogicHeightAsc);
            std::sort(decoy.fragments.rbegin(), decoy.fragments.rend(), peptideSortLogicHeightAsc);

            if (target.fragments.size() > LibraryCommon::MaxF) {
                target.fragments.resize(LibraryCommon::MaxF);
            }

            if (decoy.fragments.size() > LibraryCommon::MaxF) {
                decoy.fragments.resize(LibraryCommon::MaxF);
            }

            float max = 0.0;
            for (int i = 0; i < target.fragments.size(); i++) {
                if (target.fragments[i].height > max) {
                    max = target.fragments[i].height;
                }
            }

            for (int i = 0; i < target.fragments.size(); i++) {
                target.fragments[i].height /= max;
            }

            max = 0.0;
            for (int i = 0; i < decoy.fragments.size(); i++) {
                if (decoy.fragments[i].height > max) {
                    max = decoy.fragments[i].height;
                }
            }

            for (int i = 0; i < decoy.fragments.size(); i++) {
                decoy.fragments[i].height /= max;
            }

            target.length = decoy.length = LibraryCommon::peptide_length(name);
        }

        friend bool operator < (const Entry &left, const Entry &right) { return left.name < right.name; }

        inline void generate_decoy() {
            int i;

            decoy = target;
            auto aas = LibraryCommon::get_aas(name);
            int m = aas.size() - 2;
            if (m < 0) return;

            bool recognise = false;
            float gen_acc = 0.0;
            std::vector<Ion> pattern;
            for (auto &f : decoy.fragments) if (!(f.type & 3) || f.charge <= 0) {
                    recognise = true;
                    break;
                }
            if (recognise || LibraryCommon::ReverseDecoys) {

                pattern = recognise_fragments(
                        name,
                        target.fragments,
                        gen_acc,
                        target.charge,
                        false,
                        LibraryCommon::MaxRecCharge,
                        LibraryCommon::MaxRecLoss
                        );

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

            if (!LibraryCommon::ReverseDecoys) {
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

            auto seq = LibraryCommon::get_sequence(name, &target.no_cal);
            if (!seq.size()) return;

            float gen_acc = 0.0;

            auto pattern = recognise_fragments(
                    seq,
                    target.fragments,
                    gen_acc,
                    target.charge,
                    false,
                    LibraryCommon::MaxRecCharge,
                    LibraryCommon::MaxRecLoss
                    );

            for (i = 0; i < pattern.size(); i++) {

                target.fragments[i].charge = pattern[i].charge;

                if (i < decoy.fragments.size()) {
                    decoy.fragments[i].charge = pattern[i].charge;
                }
            }
        }

        std::vector<Ion> fragment() {
            assert(MaxF > 1000);

            int cnt = 0;
            auto seq = LibraryCommon::get_sequence(name);

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
            ReadWriteCommonFunctons::write_string(out, name);
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

        ReadWriteCommonFunctons::write_string(out, name);
        ReadWriteCommonFunctons::write_string(out, fasta_names);
        ReadWriteCommonFunctons::write_array(out, proteins);
        ReadWriteCommonFunctons::write_array(out, protein_ids);
        ReadWriteCommonFunctons::write_strings(out, precursors);
        ReadWriteCommonFunctons::write_strings(out, names);
        ReadWriteCommonFunctons::write_strings(out, genes);
        out.write((char*)&iRT_min, sizeof(double));
        out.write((char*)&iRT_max, sizeof(double));
        ReadWriteCommonFunctons::write_array(out, entries);
        ReadWriteCommonFunctons::write_vector(out, elution_groups);
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

        qDebug() << "Saving the library to " << QString::fromStdString(file_name);

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

};


#endif //PYTHIADIACPP_LIBRARYREADER_H
