//
// Created by andrewnichols on 8/8/24.
//

#include "ObjectCSVWriters.h"

Err ObjectCSVWriters::writeRawPointerToFile(
	float* vec,
	int size,
	const QString &filePathDestination
	) {

	ERR_INIT

	QVector<float> vectorFromPointerIntx(size, 0);
	std::copy_n(vec, size, vectorFromPointerIntx.data());
	e = writeVectorToFile(vectorFromPointerIntx, filePathDestination); ree;

	ERR_RETURN
}

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

Err ObjectCSVWriters::writeScanPoints(
	const QVector<PointFF> &scanPoints,
	const QString &filePathDestination
	) {

	ERR_INIT

	QFile file(filePathDestination);
	e = ErrorUtils::isTrue(file.open(QIODevice::WriteOnly)); ree;

	QTextStream stream(&file);

	for (const PointFF &v : scanPoints) {
		stream
		<< v.x() << S_GLOBAL_SETTINGS.COMMA
		<< v.y() <<  '\n';
	}

	file.close();

	ERR_RETURN
}


Err ObjectCSVWriters::readScanPoints(
	const QString &filePathSource,
	QVector<PointFF> *scanPoints
	) {

	ERR_INIT

	QFile file(filePathSource);
	e = ErrorUtils::isTrue(file.open(QIODevice::ReadOnly)); ree;

	QTextStream stream(&file);

	scanPoints->clear();

	while (!stream.atEnd()) {
		QString line = stream.readLine().trimmed();
		if (line.isEmpty()) {
			continue;
		}

		const QStringList parts = line.split(S_GLOBAL_SETTINGS.COMMA);
		e = ErrorUtils::isTrue(parts.size() == 2); ree;

		bool okX = false;
		bool okY = false;
		const float x = parts[0].toFloat(&okX);
		const float y = parts[1].toFloat(&okY);

		e = ErrorUtils::isTrue(okX && okY); ree;

		scanPoints->push_back({x, y});
	}

	file.close();

	ERR_RETURN
}