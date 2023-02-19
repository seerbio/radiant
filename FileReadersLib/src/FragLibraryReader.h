//
// Created by anichols on 2/18/23.
//

#ifndef PYTHIADIACPP_FRAGLIBRARYREADER_H
#define PYTHIADIACPP_FRAGLIBRARYREADER_H

#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReaderBase.h"

struct FILEREADERSLIB_EXPORTS FragLibIon {
    int peptideId = -1;
    double mzFrag = -1.0;
    double peptideMass = -1.0;
};

class FILEREADERSLIB_EXPORTS FragLibraryReader : public ParquetReaderBase {

public:

    FragLibraryReader() = default;
    ~FragLibraryReader() = default;

    bool readFile(const std::string &fileURI) override;

    bool writeFile(const std::string &outputFilePath) override;


private:

    std::vector<FragLibIon> m_fragLibIons;

};


#endif //PYTHIADIACPP_FRAGLIBRARYREADER_H
