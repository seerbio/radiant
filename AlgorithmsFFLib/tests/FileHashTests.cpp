#include "CandidateBundleIO.h"
#include <QCryptographicHash>
#include <QFile>
#include <QTemporaryDir>
#include <QtTest>

class FileHashTests : public QObject {
    Q_OBJECT
private slots:
    void completeByteFixtures_data();
    void completeByteFixtures();
    void knownDigests();
    void fileErrorsAndLinks();
    void virtualFile();
};

void FileHashTests::completeByteFixtures_data() {
    QTest::addColumn<int>("size");
    for (int size : {0, 1, 2, 3, 55, 56, 63, 64, 65, 127, 128, 129, 255, 256, 257,
                     1023, 1024, 1025, 4095, 4096, 4097, 65535, 65536, 65537,
                     1048575, 1048576, 1048577, 2097153}) {
        QTest::newRow(QByteArray::number(size).constData()) << size;
    }
}

void FileHashTests::completeByteFixtures() {
    QFETCH(int, size);
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    QByteArray bytes(size, Qt::Uninitialized);
    for (int i = 0; i < size; ++i) bytes[i] = char(i & 255);
    const QString expected = QString::fromLatin1(
        QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex());
    const QString path = temporary.filePath(QStringLiteral("bytes.bin"));
    QFile output(path);
    QVERIFY(output.open(QIODevice::WriteOnly));
    QCOMPARE(output.write(bytes), qint64(size));
    output.close();
    QCOMPARE(CandidateBundleIO::fileHash(path), expected);
}

void FileHashTests::knownDigests() {
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    const QString path = temporary.filePath(QString::fromUtf8("checksum-μ-质谱.bin"));
    QFile output(path);
    QVERIFY(output.open(QIODevice::WriteOnly));
    output.close();
    QCOMPARE(CandidateBundleIO::fileHash(path),
             QStringLiteral("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"));
    QVERIFY(output.open(QIODevice::WriteOnly));
    QCOMPARE(output.write("abc", 3), qint64(3));
    output.close();
    QCOMPARE(CandidateBundleIO::fileHash(path),
             QStringLiteral("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"));
}

void FileHashTests::fileErrorsAndLinks() {
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    QVERIFY(CandidateBundleIO::fileHash(temporary.filePath("missing")).isNull());
    QVERIFY(CandidateBundleIO::fileHash(temporary.path()).isNull());
#ifdef Q_OS_LINUX
    const QString target = temporary.filePath("target");
    QFile output(target);
    QVERIFY(output.open(QIODevice::WriteOnly));
    QCOMPARE(output.write("abc", 3), qint64(3));
    output.close();
    const QString link = temporary.filePath("link");
    QVERIFY(QFile::link(target, link));
    QCOMPARE(CandidateBundleIO::fileHash(link), CandidateBundleIO::fileHash(target));
    QVERIFY(QFile::remove(target));
    QVERIFY(CandidateBundleIO::fileHash(link).isNull());
#endif
}

void FileHashTests::virtualFile() {
#ifdef Q_OS_LINUX
    // Its reported size is zero; hashing must still read every available byte.
    QFile file(QStringLiteral("/proc/version"));
    QVERIFY(file.open(QIODevice::ReadOnly));
    const QByteArray bytes = file.readAll();
    QVERIFY(!bytes.isEmpty());
    QCOMPARE(CandidateBundleIO::fileHash(file.fileName()),
             QString::fromLatin1(QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex()));
#endif
}

QTEST_MAIN(FileHashTests)
#include "FileHashTests.moc"
