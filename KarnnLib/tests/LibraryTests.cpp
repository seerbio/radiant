//
// Created by anichols on 11/07/2021.
//

#include "Library.h"

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

    std::stringstream ls(diannLibraryFile.toStdString());
    Library lib;
    lib.load(diannLibraryFile.toStdString().c_str());

    for (const Library::Entry &entry : lib.entries) {

        qDebug() << QString::fromStdString(entry.name) << entry.target.iRT << entry.target.sRT;
        for (const Product &pr : entry.target.fragments) {
            qDebug() << pr.mz << pr.height << toascii(pr.charge) << toascii(pr.index)
                     << toascii(pr.type) << QString::fromStdString(getLoss(static_cast<int>(toascii(pr.loss))));
        }

    }

    qDebug() << lib.entries.size();


}




QTEST_MAIN(LibraryTests)
#include "LibraryTests.moc"
