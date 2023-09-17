//
// Created by anichols on 11/07/2021.
//

#include "LibraryReaderr.h"

#include <QDebug>
#include <QtTest/QtTest>


class LibraryTestss : public QObject
{
    Q_OBJECT

public:
    LibraryTestss() = default;
    ~LibraryTestss() override = default;

private Q_SLOTS:

    void readLibrary();

};

void LibraryTestss::readLibrary() {

    LibraryReaderr r;
    r.test();

}


QTEST_MAIN(LibraryTestss)
#include "LibraryTestss.moc"
