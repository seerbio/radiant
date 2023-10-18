#include "BoostRTreeWrapper.h"

#include <QtTest/QtTest>


class BoostRTreeWrapperTests : public QObject
{
    Q_OBJECT

public:
    BoostRTreeWrapperTests() = default;
    ~BoostRTreeWrapperTests() override = default;

private Q_SLOTS:
    void initPointTree();

};

void BoostRTreeWrapperTests::initPointTree() {

    ERR_INIT

    const QVector<RTreePointData2D> points = {
            RTreePointData2D(1.0, 1.0, 1.0),
            RTreePointData2D(2.0, 2.0, 2.0),
            RTreePointData2D(3.0, 3.0, 3.0),
            RTreePointData2D(4.0, 4.0, 4.0)
    };

    BoostRTreeWrapper boostRTreeWrapper;
    e = boostRTreeWrapper.init(points);
    QCOMPARE(e, eNoError);

    QVector<RTreePointData2D> returnedPoints;
    e = boostRTreeWrapper.getPoints(1.5, 3.5, 1.5, 3.5, &returnedPoints);

    const RTreePointData2D &p1 = returnedPoints.front();
    QCOMPARE(p1.x, 2.0);
    QCOMPARE(p1.y, 2.0);
    QCOMPARE(p1.val, 2.0);
    
    const RTreePointData2D &p2 = returnedPoints.back();
    QCOMPARE(p2.x, 3.0);
    QCOMPARE(p2.y, 3.0);
    QCOMPARE(p2.val, 3.0);

}


QTEST_MAIN(BoostRTreeWrapperTests)
#include "BoostRTreeWrapperTests.moc"
