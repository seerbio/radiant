//
// Created by anichols on 11/07/2021.
//

#include "LibraryReader.h"

#include "LibraryCommon.h"

#include <QDebug>
#include <QtTest/QtTest>


class LibraryTests : public QObject
{
    Q_OBJECT

public:
    LibraryTests() = default;
    ~LibraryTests() override = default;


private Q_SLOTS:

    void readLibrary();


};


void LibraryTests::readLibrary() {

    const QString diannLibraryFile = "/home/anichols/Desktop/Testing/lib.predicted.speclib";

    LibraryReader lib;
    lib.load(diannLibraryFile.toStdString().c_str());

    for (const Entry &entry : lib.getEntries()) {

        qDebug() << QString::fromStdString(entry.name) << entry.target.iRT << entry.target.sRT;
        for (const Product &pr : entry.target.fragments) {
            qDebug() << pr.mz << pr.height << toascii(pr.charge) << toascii(pr.index)
                     << toascii(pr.type) << QString::fromStdString(LibraryCommon::getLoss(static_cast<int>(toascii(pr.loss))));
        }

    }


}




QTEST_MAIN(LibraryTests)
#include "LibraryTests.moc"
