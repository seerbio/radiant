//
// Created by anichols on 11/07/2021.
//

#include "MsReaderBase.h"

#include <QtTest/QtTest>


class MsReaderBaseTests : public QObject
{
    Q_OBJECT

public:
    MsReaderBaseTests() = default;
    ~MsReaderBaseTests() override = default;

private Q_SLOTS:





private:
    //TODO use proper path procedures.
    const QString m_filepath
            = QStringLiteral("/home/anichols/Desktop/RawData/EXP21155_2021ms0609X7_A.raw.mzML");

};


QTEST_MAIN(MsReaderBaseTests)
#include "MsReaderBaseTests.moc"
