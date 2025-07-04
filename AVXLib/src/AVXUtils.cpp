//
// Created by andrewnichols on 7/4/25.
//

#include "AVXUtils.h"


size_t AVXUtils::calculateNextAlignedBlockSize(
	size_t arrSize,
	size_t byteSize
	) {

	while (arrSize % byteSize != 0) {
		arrSize++;
	}

	return arrSize;
}
