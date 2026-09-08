//
// Created by Codex on 9/8/26.
//

#ifndef IMHANDLINGMODE_H
#define IMHANDLINGMODE_H

#include "FileReadersLib_Exports.h"

#include <QString>

enum class ImHandlingMode {
    Centroid,
    Summed,
    Raw4D,
};

FILEREADERSLIB_EXPORTS QString imHandlingModeToString(ImHandlingMode mode);
FILEREADERSLIB_EXPORTS bool imHandlingModeFromString(
    const QString &text,
    ImHandlingMode *mode
    );

#endif // IMHANDLINGMODE_H
