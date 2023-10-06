//
// Created by anichols on 7/22/23.
//

#ifndef PYTHIADIACPP_IRTPREDICTRON_H
#define PYTHIADIACPP_IRTPREDICTRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;


class ALGORITHMSLIB_EXPORTS IRTPredictron {

public:

    IRTPredictron();
    ~IRTPredictron();

    Err init(const QString &modelFilePath);

    Err batchPredictIRT(
            const QStringList &peptideStringWithModsList,
            QVector<float> *rawPredictionResults
    );

private:

    Q_DISABLE_COPY(IRTPredictron) class Private;
    const QScopedPointer<Private> d_ptr;



};


#endif //PYTHIADIACPP_IRTPREDICTRON_H
