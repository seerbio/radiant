//
// Created by andrewnichols on 8/15/24.
//

#ifndef QUANFILEBUILDER_H
#define QUANFILEBUILDER_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"

class CandidateScores;
class MsFrame;

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS QuanFileBuilder {

public:

    static Err buildQuanFile(
        const QVector<CandidateScores*> &candidateScores,
        const QString& outputFilePath
        );



};



#endif //QUANFILEBUILDER_H
