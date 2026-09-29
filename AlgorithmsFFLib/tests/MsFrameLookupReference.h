#pragma once

#include "MsFrame.h"
#include <Eigen/Dense>
#include <nanoflann.hpp>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

// Independent map and nearest-neighbor behavior from MsFrame at commit
// 2b6001a874ef09e7a884975e0e7d538e733fa6b9. This retains the original dynamic
// query buffers and map lookups; it does not call the production cache.
class MsFrameLookupReference {
    using KDTree = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd>;
    QMap<FrameIndex, ScanNumber> m_frames;
    QMap<ScanNumber, ScanTime> m_times;
    QMap<Index, FrameIndex> m_indexToFrame;
    Eigen::MatrixXd m_matrix;
    std::unique_ptr<KDTree> m_tree;

public:
    void init(const QList<ScanNumber> &numbers, const QMap<ScanNumber, ScanTime> &times) {
        m_times = times;
        for (const ScanNumber number : numbers)
            m_frames.insert(m_frames.size(), number);
        m_tree.reset();
        m_matrix.resize(m_frames.size(), 2);
        int index = 0;
        for (auto it = m_frames.begin(); it != m_frames.end(); ++it) {
            m_matrix.coeffRef(index, 0) = m_times.value(it.value());
            m_matrix.coeffRef(index, 1) = 0.0;
            m_indexToFrame.insert(index++, it.key());
        }
        m_tree = std::make_unique<KDTree>(2, m_matrix, 15);
    }

    int scanCount() const { return m_frames.size(); }

    ScanNumber scanNumberFromFrameIndex(FrameIndex frame) const {
        if (frame > m_frames.lastKey()) return m_frames.last();
        if (frame < 1) return m_frames.first();
        return m_frames.value(frame);
    }

    ScanTime scanTimeFromScanNumber(ScanNumber number) const {
        return m_times.value(number);
    }

    ScanTime scanTimeFromFrameIndex(FrameIndex frame) const {
        return scanTimeFromScanNumber(scanNumberFromFrameIndex(frame));
    }

    FrameIndex frameIndexFromScanTime(ScanTime time) const {
        const size_t numResults = 1;
        std::vector<double> query = {time, 0.0};
        std::vector<long> result(numResults);
        std::vector<double> distance(numResults);
        m_tree->index->knnSearch(query.data(), numResults, result.data(), distance.data());
        return m_indexToFrame.value(static_cast<int>(result.front()));
    }
};

struct FrameLookupInputs {
    QMap<ScanNumber, ScanPoints> points;
    QMap<ScanNumber, ScanTime> times;
};

inline FrameLookupInputs makeFrameLookupInputs(int count, int kind) {
    FrameLookupInputs inputs;
    for (int i = 0; i < count; ++i) {
        const int number = kind == 4 ? -1000000000 + i * (2000000000 / count)
                         : kind == 5 ? std::numeric_limits<int>::max() - count * 3 + i * 3
                         : kind == 6 ? std::numeric_limits<int>::min() + i * 3
                         : kind == 7 ? 16777200 + i * 3 : i * 21 - 25;
        float time = float(i) / 60.f;
        if (kind == 1) time = float(i / 3) / 60.f;
        if (kind == 2) time = float((i * 17) % count - count / 2) / 60.f;
        if (kind == 3) time = i % 3 == 0 ? -0.f : float(i % 3) * std::numeric_limits<float>::denorm_min();
        inputs.points.insert(number, {{100.f, 1.f}});
        if (kind != 8 || i % 3 != 0) inputs.times.insert(number, time);
    }
    inputs.times.insert(inputs.points.lastKey() + 1, -0.f);
    return inputs;
}

inline std::uint32_t scanTimeBits(ScanTime time) {
    static_assert(sizeof(time) == sizeof(std::uint32_t));
    std::uint32_t bits;
    std::memcpy(&bits, &time, sizeof(bits));
    return bits;
}
