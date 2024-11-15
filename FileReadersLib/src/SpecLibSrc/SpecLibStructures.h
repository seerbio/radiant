//
// Created by andrewnichols on 11/15/24.
//

#ifndef SPECLIB_STRUCTS_H
#define SPECLIB_STRUCTS_H

#include "FileReadersLib_Exports.h"
#include "Readers.h"

#include <set>
#include <vector>
#include <string>
#include <regex>


class FILEREADERSLIB_EXPORTS Isoform {

public:
	std::string id;
	mutable std::string name;
	mutable std::string gene;
	mutable std::set<int> precursorIndexes;
	int nameIndex = 0;
	int geneIndex = 0;
	bool swissprot = true;

	template <class F> void
	read(F &input) {
		int sp = 0;
		int size = 0;
		input.read(reinterpret_cast<char*>(&sp), sizeof(int));
		input.read(reinterpret_cast<char*>(&size), sizeof(int));
		swissprot = sp;

		Readers::readString(input, id);
		Readers::readString(input, name);
		Readers::readString(input, gene);
		input.read(reinterpret_cast<char*>(&nameIndex), sizeof(int));
		input.read(reinterpret_cast<char*>(&geneIndex), sizeof(int));
		precursorIndexes.clear();
		for (int i = 0; i < size; i++) {
			int pr = -1;
			input.read(reinterpret_cast<char*>(&pr), sizeof(int));
			if (pr >= 0) precursorIndexes.insert(pr);
		}
	}
};

class FILEREADERSLIB_EXPORTS PG {

public:
	std::string ids;
	mutable std::string names;
	mutable std::string genes;
	mutable std::vector<int> precursorIndexes; // precursor indices in the library
	mutable std::set<int> proteins;
	std::vector<int> nameIndices;
	std::vector<int> geneIndices;

	template <class F> void
	read(F &input) {
		int sp = 0;
		int sizeP = 0;
		input.read(reinterpret_cast<char*>(&sizeP), sizeof(int));

		Readers::readString(input, ids);
		Readers::readString(input, names);
		Readers::readString(input, genes);
		Readers::readVector(input, nameIndices);
		Readers::readVector(input, geneIndices);
		Readers::readVector(input, precursorIndexes);
		for (int i = 0; i < sizeP; i++) {
			int p = -1;
			input.read(reinterpret_cast<char*>(&p), sizeof(int));
			if (p >= 0) proteins.insert(p);
		}
	}
};

class FILEREADERSLIB_EXPORTS Product {
public:
	float mz = 0.0;
	float height = 0.0;
	uint8_t charge = 0;
	uint8_t type = 0;
	uint8_t index = 0;
	uint8_t loss = 0;
};

class FILEREADERSLIB_EXPORTS Peptide {
public:
	int index = 0;
	int charge = 0;
	int length = 0;

	float mz = 0.0;
	float iRT = 0.0;
	float sRT = 0.0;
	float iIM = 0.0;
	float sIM = 0.0;
	float libQvalue = 0.0;

	std::vector<Product> fragments;

	template <class F>
	void read(F &input, int version) {
		input.read(reinterpret_cast<char*>(&index), sizeof(int));
		input.read(reinterpret_cast<char*>(&charge), sizeof(int));
		input.read(reinterpret_cast<char*>(&length), sizeof(int));

		input.read(reinterpret_cast<char*>(&mz), sizeof(float));
		input.read(reinterpret_cast<char*>(&iRT), sizeof(float));
		input.read(reinterpret_cast<char*>(&sRT), sizeof(float));

		if (version <= -2) {
			input.read(reinterpret_cast<char*>(&libQvalue), sizeof(float));
			input.read(reinterpret_cast<char*>(&iIM), sizeof(float));
			input.read(reinterpret_cast<char*>(&sIM), sizeof(float));
		}

		Readers::readVector(input, fragments);
	}
};

#endif //SPECLIB_STRUCTS_H
