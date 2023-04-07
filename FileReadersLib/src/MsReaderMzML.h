#ifndef MSREADERMZML_H
#define MSREADERMZML_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "MsReaderBase.h"

#include <QByteArray>
#include <QPointF>
#include <QScopedPointer>
#include <QSharedPointer>
#include <QVector>
#include <QXmlStreamReader>


using namespace Error;


class FILEREADERSLIB_EXPORTS MsReaderMzML : public MsReaderBase {

public:

    MsReaderMzML();

    ~MsReaderMzML();

    Err openFile(const QString &filePath) override;

    Err closeFile() override;

    MsReaderBase msReaderBase();

private:

    Q_DISABLE_COPY(MsReaderMzML)
    class PrivateData;
    QScopedPointer<PrivateData> m_d;
};


#endif // MSREADERMZML_H
