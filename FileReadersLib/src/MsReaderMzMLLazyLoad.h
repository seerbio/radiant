#ifndef MSREADERMZMLLAZYLOAD_H
#define MSREADERMZMLLAZYLOAD_H

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

namespace MzMLLazyLoadNamespace {
    const QString MZML_SUFFIX = QStringLiteral("mzML");
}

class FILEREADERSLIB_EXPORTS MsReaderMzMLLazyLoad : public MsReaderBase {

public:

    friend class MsReaderMZMLLazyLoadTests;
    friend class TurboXICTests;

    MsReaderMzMLLazyLoad();

    ~MsReaderMzMLLazyLoad() override;

    /**
    * @brief Opens an MzML file in MsReaderMzMLMapped::PrivateData.
    *
    * This function initializes the opening of an MzML file specified by the `filename`
    * parameter in the PrivateData of MsReaderMzMLMapped. It first attempts to close any existing
    * file using the `closeFile` method. Then, it creates a new QFile instance and opens
    * the file for reading. Finally, it initializes a QXmlStreamReader to parse the XML content
    * of the file using the `parse` method. Any encountered errors during the process are
    * indicated by the returned Err code.
    *
    * @param filename The name of the MzML file to be opened.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err openFile(const QString &filePath) override;

    /**
    * @brief Closes the MzML file in MsReaderMzMLMapped::PrivateData.
    *
    * This function clears the member variables `m_msScanInfo` and `m_scanPoints`,
    * ensuring that they are empty, and resets the QFile instance to close the file.
    * Any encountered errors during the process are indicated by the returned Err code.
    *
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err closeFile() override;

    Err getMzTargetScanPoints(
        const MzTargetKey& targetKey,
        QVector<MsScan> *msScans
        ) override;

    /**
    * @brief This function sets scan information and scan points of an MsReaderBase instance from an MsReaderMzMLMapped instance.
    *
    * First, it initializes an MsReaderBase instance. Then, it retrieves the scan information and scan points from
    * the current MsReaderMzMLMapped instance and sets them to the newly created MsReaderBase instance. The function
    * highlights the fact that setMsScanInfo should be called before setScanPoints to avoid runtime errors.
    * Lastly, the function returns the newly created MsReaderBase instance.
    *
    * @return MsReaderBase - an instance of MsReaderBase with scan information and scan points set from the current
    * MsReaderMzMLMapped instance.
    */
    // MsReaderBase msReaderBase() const;

    Err extractScanPoints(
	    const QVector<MsScanInfo*> &msScanInfos,
		QVector<MsScan> *msScans
        ) override;

	static Err extractScanPoints(
		const QString &fileName,
		const QVector<MsScanInfo*> &msScanInfos,
		QVector<MsScan> *msScans
		);

	// static Err extractScanPoints(
	// 	const QString &fileName,
	// 	QVector<MsScanInfo> *msScanInfos,
	// 	QVector<MsScan> *msScans
	// 	);

private:

    Q_DISABLE_COPY(MsReaderMzMLLazyLoad)
    class PrivateData;
    QScopedPointer<PrivateData> m_d;

};


#endif // MSREADERMZMLLAZYLOAD_H
