#ifndef MSREADERTHERMORAW_H
#define MSREADERTHERMORAW_H

#include "FileReadersLib_Exports.h"
#include "MsReaderBase.h"

/**
 * Reads Thermo RAW spectra through the RadiantThermoRawReader binary bridge.
 *
 * The bridge uses the ThermoRawFileParser and Thermo RawFileReader assemblies
 * directly. Spectra and metadata are streamed as binary data; no mzML is made.
 */
class FILEREADERSLIB_EXPORTS MsReaderThermoRaw : public MsReaderBase {
public:
    explicit MsReaderThermoRaw(const QString &readerPath = QString());
    ~MsReaderThermoRaw() override = default;

    using MsReaderBase::openFile;
    Err openFile(const QString &filePath) override;
    Err closeFile() override;

private:
    Q_DISABLE_COPY(MsReaderThermoRaw)

    QString m_readerPath;
};

#endif // MSREADERTHERMORAW_H
