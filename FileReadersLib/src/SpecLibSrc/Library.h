//
// Created by andrewnichols on 11/15/24.
//

#ifndef LIBRARY_H
#define LIBRARY_H

#include "ChemConstants.h"
#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "FragLibReaderRow.h"
#include "PeptideStringWithMods.h"
#include "Readers.h"
#include "SpecLibStructures.h"

#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>


using namespace Error;


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

	static Err extractIonLabel(
		unsigned int ionIndexRaw,
		unsigned int ionTypeInt,
		unsigned int ionCharge,
		unsigned int ionLoss,
		int peptideLen,
		IonLabel *ionLabel
		) {

		ERR_INIT

		e = ErrorUtils::isWithinRange(
			static_cast<int>(ionTypeInt),
			1, 2,
			ErrorUtilsParam::IncludeThreshold,
			eValueError
			); ree;

		e = ErrorUtils::isEqual(static_cast<int>(ionLoss), 0); ree;

		const QChar ionType = ionTypeInt == 1 ? 'b' : 'y';
		const int ionIndex = ionTypeInt == 1 ? ionIndexRaw : peptideLen - ionIndexRaw;

		IonLabel ionLabelLocal = ionType + QString::number(ionIndex);
		ionLabelLocal = ionCharge > 1 ? ionLabelLocal + '^' + QString::number(static_cast<int>(ionCharge)) : ionLabelLocal;

		*ionLabel = ionLabelLocal;

		ERR_RETURN;
	}

	static QPair<Err, QVector<FragLibReaderRow>> entriesLoadLogic(
		const QVector<int> &indexes,
		const std::vector<Entry> &entries,
		const std::vector<std::string> &precursors,
		const std::vector<PG> &proteinIds
		) {

		ERR_INIT

		QVector<FragLibReaderRow> fragLibReaderRows;
		fragLibReaderRows.reserve(indexes.size());

		for (int index : indexes) {
			const Entry &entry = entries[index];

			const auto &[target, decoy, entryFlags, proteotypic, name, pidIndex, pgQValue, ptmQvalue, siteConf] = entry;
			const auto &p = precursors.at(index);

			FragLibReaderRow flrr;

			const QString proteinGroups = QString::fromStdString(proteinIds.at(entry.pidIndex).names);
			flrr.proteinGroups = proteinGroups;

			QString specLibPeptide = QString::fromStdString(p);
			const QChar specLibPeptideCharge = specLibPeptide.back();

			specLibPeptide.chop(1);
			const PeptideString peptideString = PeptideStringWithMods(specLibPeptide).removeUniModChars();
			specLibPeptide += '|' + specLibPeptideCharge;

			flrr.peptideSequenceChargeKey = specLibPeptide;
			e = ErrorUtils::toInt(QString(specLibPeptideCharge), &flrr.precursorCharge); rree;
			flrr.iM = target.iIM;
			flrr.iRT = target.iRT;
			flrr.isDecoy = false;

			const int fragmentsSize = static_cast<int>(target.fragments.size());
			flrr.intensityVals.reserve(fragmentsSize);
			flrr.mzVals.reserve(fragmentsSize);
			flrr.mass = (target.mz * static_cast<float>(flrr.precursorCharge)) - (ChemConstants::PROTON * flrr.precursorCharge);

			QStringList ionLabelsList;
			ionLabelsList.reserve(fragmentsSize);

			for (int i = 0; i < fragmentsSize; i++) {
				const auto ionIndexRaw = static_cast<unsigned int>(target.fragments.at(i).index);
				const auto ionTypeInt = static_cast<unsigned int>(target.fragments.at(i).type);
				const auto ionCharge = static_cast<unsigned int>(target.fragments.at(i).charge);
				const auto ionLoss = static_cast<unsigned int>(target.fragments.at(i).loss);
				const float intensityVal = target.fragments.at(i).height;
				const float mzVal = target.fragments.at(i).mz;

				if (
					constexpr float intensityCutoff = 0.01;
					ionCharge > flrr.precursorCharge || intensityVal < intensityCutoff)
					{
					continue;
				}

				flrr.mzVals.push_back(mzVal);
				flrr.intensityVals.push_back(intensityVal);

				IonLabel ionLabel;
				e = extractIonLabel(
					ionIndexRaw,
					ionTypeInt,
					ionCharge,
					ionLoss,
					peptideString.size(),
					&ionLabel
					); rree;
				ionLabelsList.push_back(ionLabel);
			}

			flrr.ionLabels = ionLabelsList.join(S_GLOBAL_SETTINGS.SEPARATOR);
			fragLibReaderRows.push_back(flrr);
		}

		return {e, fragLibReaderRows};
	}

	template<class F>
	Err read(
		F &input,
		QVector<FragLibReaderRow> *fragLibReaderRows
		) {

		ERR_INIT

		int gd = 0;
		int gc = 0;
		int ip = 0;
		int version = 0;

		input.read(reinterpret_cast<char*>(&version), sizeof(int));
		if (version >= 0) gd = version;
		else input.read(reinterpret_cast<char*>(&gd), sizeof(int));

		if (version < -3) {
			std::cout << "ERROR: version is not supported" << -version << std::endl;
			rrr(eFileError);
		}

		input.read(reinterpret_cast<char*>(&gc), sizeof(int));
		input.read(reinterpret_cast<char*>(&ip), sizeof(int));
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

		genDecoys = gd;
		genCharges = gc;
		inferProteotypicity = ip;

		if (version <= -1 && input.peek() != std::char_traits<char>::eof()) {
			Readers::readVector(input, elution_groups);
		}

		e = ErrorUtils::isNotEmpty(entries); ree;
		fragLibReaderRows->reserve(static_cast<int>(entries.size()));

#define ENTRIES_PARALLEL
#ifdef ENTRIES_PARALLEL
		QVector<int> indexes(entries.size());
		std::iota(indexes.begin(), indexes.end(), 0);

		QVector<QVector<int>> indexesTranched;
		e = ParallelUtils::trancheVectorForParallelizationInOrder(
			indexes,
			ParallelUtils::numberOfAvailableSystemProcessors(),
			0,
			&indexesTranched
			); ree;

		const auto binderLogic = std::bind(
			entriesLoadLogic,
			std::placeholders::_1,
			entries,
			precursors,
			proteinIds
			);

		QFuture<QPair<Err, QVector<FragLibReaderRow>>> futures = QtConcurrent::mapped(
			indexesTranched,
			binderLogic
			);
		futures.waitForFinished();

		for (const QPair<Err, QVector<FragLibReaderRow>> &pr : futures) {
			e = pr.first; ree;
			fragLibReaderRows->append(pr.second);
		}
#else
		for (int i = 0; i < entries.size(); i++) {
			QPair<Err, FragLibReaderRow> pr = entriesLoadLogic(
				i,
				entries,
				precursors,
				proteinIds
				);

			e = pr.first; ree;
			fragLibReaderRows->push_back(pr.second);
		}
#endif
		ERR_RETURN
	}
};


#endif //LIBRARY_H
