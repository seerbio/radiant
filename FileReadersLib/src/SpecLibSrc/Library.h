//
// Created by andrewnichols on 11/15/24.
//

#ifndef LIBRARY_H
#define LIBRARY_H

#include "FileReadersLib_Exports.h"
#include "Readers.h"
#include "SpecLibStructures.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <regex>

class FILEREADERSLIB_EXPORTS Entry {

public:
	Peptide target;
	Peptide decoy;
	int entryFlags = 0;
	int proteotypic = 0;
	std::string name;
	int pidIndex = 0;
	float pgQValue = 0.0;
	float ptmQvalue = 1.0;
	float siteConf = 1.0;

	template <class F>
	void read(F &input, int v) {

		target.read(input, v);

		int dc = 0;
		input.read(reinterpret_cast<char*>(&dc), sizeof(int));

		if (dc) {
			decoy.read(input, v);
		}

		int ff = 0;
		int prt = 0;
		input.read(reinterpret_cast<char*>(&ff), sizeof(int));
		input.read(reinterpret_cast<char*>(&prt), sizeof(int));
		entryFlags = ff, proteotypic = prt;
		input.read(reinterpret_cast<char*>(&pidIndex), sizeof(int));
		Readers::readString(input, name);
		if (v <= -3) {
			input.read(reinterpret_cast<char*>(&pgQValue), sizeof(float));
			input.read(reinterpret_cast<char*>(&ptmQvalue), sizeof(float));
			input.read(reinterpret_cast<char*>(&siteConf), sizeof(float));
		}
	}
};


class FILEREADERSLIB_EXPORTS Library {

public:
	std::string name, fastaNames;
	std::vector<Isoform> proteins;
	std::vector<PG> proteinIds;
	std::vector<PG> proteinGroups;
	std::vector<PG> geneGroups;
	std::vector<int> ggIndex;
	std::vector<std::string> precursors;
	std::vector<std::string> names;
	std::vector<std::string> genes;

	std::vector<int> elution_groups;

	double iRTMin = 0.0;
	double iRTMax = 0.0;
	bool genDecoys = true;
	bool genCharges = true;
	bool inferProteotypicity = true;
	bool fromSpeclib = false;

	std::vector<Entry> entries;

	template<class F> void
	read(F &input) {
		int gd = 0;
		int gc = 0;
		int ip = 0;
		int version = 0;

		input.read(reinterpret_cast<char*>(&version), sizeof(int));
		if (version >= 0) gd = version;
		else input.read(reinterpret_cast<char*>(&gd), sizeof(int));

		if (version < -3) {
			std::cout << "ERROR: version " << -version << " of the .speclib format is not supported by this DIA-NN version" << std::endl;
			exit(-1);
		}

		input.read(reinterpret_cast<char*>(&gc), sizeof(int));
		input.read(reinterpret_cast<char*>(&ip), sizeof(int));
		genDecoys = gd;
		genCharges = gc;
		inferProteotypicity = ip;

		Readers::readString(input, name);
		Readers::readString(input, fastaNames);
		Readers::readArray(input, proteins);
		Readers::readArray(input, proteinIds);
		Readers::readStrings(input, precursors);
		Readers::readStrings(input, names);
		Readers::readStrings(input, genes);
		input.read(reinterpret_cast<char*>(&iRTMin), sizeof(double));
		input.read(reinterpret_cast<char*>(&iRTMax), sizeof(double));
		Readers::readArray(input, entries, version);

		if (version <= -1 && input.peek() != std::char_traits<char>::eof()) {
			Readers::readVector(input, elution_groups);
		}

		for (int i = 0; i < entries.size(); i++) {
			const auto &[target, decoy, entryFlags, proteotypic, name, pidIndex, pgQValue, ptmQvalue, siteConf] = entries.at(i);
			const auto &p = precursors.at(i);
			std::cout
			<< p
			<< " " << target.iRT
			<< " " << target.iIM
			<< " " << static_cast<unsigned int>(target.fragments.at(0).index)
			<< " " << static_cast<unsigned int>(target.fragments.at(0).type)
			<< " " << static_cast<unsigned int>(target.fragments.at(0).charge)
			<< " " << static_cast<unsigned int>(target.fragments.at(0).loss)
			<< " " << target.fragments.at(0).height
			<< " " << target.fragments.at(0).mz
			<< std::endl;;
		}

	}
};


#endif //LIBRARY_H
