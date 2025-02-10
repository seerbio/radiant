//
// Created by andrewnichols on 8/8/24.
//

#include "ObjectCSVWriters.h"

Err ObjectCSVWriters::writeScanPoints(
    const QMap<int, QVector<PointFF>> &scanPoints,
    const QString &filePathDestination
    ) {

    ERR_INIT

    QFile file(filePathDestination);
    e = ErrorUtils::isTrue(file.open(QIODevice::WriteOnly)); ree;

    QTextStream stream(&file);

    for (auto it = scanPoints.begin(); it != scanPoints.end(); ++it) {
        const int scanNumber = it.key();
        const QVector<PointFF> &scanPointsVec = it.value();

        for (const PointFF &v : scanPointsVec) {
            stream
            << scanNumber << S_GLOBAL_SETTINGS.COMMA
            << v.x() << S_GLOBAL_SETTINGS.COMMA
            << v.y() <<  '\n';
        }

    }

    file.close();

    ERR_RETURN

}
