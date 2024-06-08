//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsReaderBase.h"

#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>

#include "GlobalSettings.h"
#include "MsReaderBase.h"


class PlayGroundTests : public QObject
{
    Q_OBJECT

public:
    PlayGroundTests() = default;
    ~PlayGroundTests() override = default;

private Q_SLOTS:

void testme();

};

void PlayGroundTests::testme() {



}

QTEST_MAIN(PlayGroundTests)
#include "PlayGroundTests.moc"
