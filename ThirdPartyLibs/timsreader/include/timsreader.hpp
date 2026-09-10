#ifndef TIMSREADER_HPP
#define TIMSREADER_HPP

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "timsreader.h"

namespace timsreader {

inline void throw_if_error(tr_status_t status) {
  if (status == TR_STATUS_OK || status == TR_STATUS_END) {
    return;
  }
  throw std::runtime_error(tr_last_error_message());
}

enum class Representation {
  DirectSubscans = TR_DIRECT_SUBSCANS,
  FlattenedWindowSpectrum = TR_FLATTENED_WINDOW_SPECTRUM,
  CentroidedFrameSpectrum = TR_CENTROIDED_FRAME_SPECTRUM,
};

struct Selector {
  tr_selector_t raw{};

  static Selector ms1() {
    Selector out;
    out.raw.kind = TR_SELECTOR_MS1;
    out.raw.isolation_window_id = 0;
    return out;
  }

  static Selector ms2() {
    Selector out;
    out.raw.kind = TR_SELECTOR_MS2;
    out.raw.isolation_window_id = 0;
    return out;
  }

  // Direct-mode helper for one logical MS2 window. Processed MS2 plans use
  // ms2() and route rows by the exported isolation_window_id column.
  static Selector ms2_window(uint32_t isolation_window_id) {
    Selector out;
    out.raw.kind = TR_SELECTOR_MS2_WINDOW;
    out.raw.isolation_window_id = isolation_window_id;
    return out;
  }
};

struct IsolationWindowInfo {
  uint32_t isolation_window_id;
  int32_t window_group_id;
  uint16_t scan_start_index;
  uint16_t scan_end_index_exclusive;
  float isolation_target_mz;
  float isolation_lower_offset;
  float isolation_upper_offset;
};

struct BorrowedBatchView {
  const ArrowArray* array = nullptr;
  const ArrowSchema* schema = nullptr;
};

inline tr_stream_config_t default_stream_config(Selector selector) {
  return tr_stream_config_default_for_selector(selector.raw);
}

class OwnedBatch {
 public:
  OwnedBatch() = default;
  explicit OwnedBatch(tr_owned_batch_t* handle) : handle_(handle) {}

  OwnedBatch(OwnedBatch&& other) noexcept : handle_(std::exchange(other.handle_, nullptr)) {}
  OwnedBatch& operator=(OwnedBatch&& other) noexcept {
    if (this != &other) {
      reset();
      handle_ = std::exchange(other.handle_, nullptr);
    }
    return *this;
  }

  OwnedBatch(const OwnedBatch&) = delete;
  OwnedBatch& operator=(const OwnedBatch&) = delete;

  ~OwnedBatch() { reset(); }

  explicit operator bool() const { return handle_ != nullptr; }

  BorrowedBatchView export_view() const {
    const ArrowArray* array = nullptr;
    const ArrowSchema* schema = nullptr;
    throw_if_error(tr_owned_batch_export(handle_, &array, &schema));
    return BorrowedBatchView{array, schema};
  }

 private:
  void reset() {
    if (handle_ != nullptr) {
      tr_owned_batch_destroy(handle_);
      handle_ = nullptr;
    }
  }

  tr_owned_batch_t* handle_ = nullptr;
};

class Stream {
 public:
  Stream() = default;
  explicit Stream(tr_stream_t* handle) : handle_(handle) {}

  Stream(Stream&& other) noexcept : handle_(std::exchange(other.handle_, nullptr)) {}
  Stream& operator=(Stream&& other) noexcept {
    if (this != &other) {
      reset();
      handle_ = std::exchange(other.handle_, nullptr);
    }
    return *this;
  }

  Stream(const Stream&) = delete;
  Stream& operator=(const Stream&) = delete;

  ~Stream() { reset(); }

  std::optional<BorrowedBatchView> next_borrowed() {
    const ArrowArray* array = nullptr;
    const ArrowSchema* schema = nullptr;
    const tr_status_t status = tr_stream_next_borrowed(handle_, &array, &schema);
    if (status == TR_STATUS_END) {
      return std::nullopt;
    }
    throw_if_error(status);
    return BorrowedBatchView{array, schema};
  }

  std::optional<OwnedBatch> next_owned() {
    tr_owned_batch_t* batch = nullptr;
    const tr_status_t status = tr_stream_next_owned(handle_, &batch);
    if (status == TR_STATUS_END) {
      return std::nullopt;
    }
    throw_if_error(status);
    return OwnedBatch(batch);
  }

 private:
  void reset() {
    if (handle_ != nullptr) {
      tr_stream_destroy(handle_);
      handle_ = nullptr;
    }
  }

  tr_stream_t* handle_ = nullptr;
};

class Plan {
 public:
  Plan() = default;
  explicit Plan(tr_plan_t* handle) : handle_(handle) {}

  Plan(Plan&& other) noexcept : handle_(std::exchange(other.handle_, nullptr)) {}
  Plan& operator=(Plan&& other) noexcept {
    if (this != &other) {
      reset();
      handle_ = std::exchange(other.handle_, nullptr);
    }
    return *this;
  }

  Plan(const Plan&) = delete;
  Plan& operator=(const Plan&) = delete;

  ~Plan() { reset(); }

  size_t row_count() const {
    size_t out = 0;
    throw_if_error(tr_plan_row_count(handle_, &out));
    return out;
  }

  tr_stream_config_t default_stream_config() const {
    tr_stream_config_t config{};
    throw_if_error(tr_plan_default_stream_config(handle_, &config));
    return config;
  }

  // Window ids that may appear in this partition's Arrow rows.
  std::vector<uint32_t> isolation_window_ids() const {
    size_t count = 0;
    throw_if_error(tr_plan_isolation_window_count(handle_, &count));
    std::vector<uint32_t> out;
    out.reserve(count);
    for (size_t index = 0; index < count; ++index) {
      uint32_t value = 0;
      throw_if_error(tr_plan_get_isolation_window_id(handle_, index, &value));
      out.push_back(value);
    }
    return out;
  }

  std::vector<Plan> split(size_t count) const {
    size_t actual = 0;
    throw_if_error(tr_plan_partition_count(handle_, count, &actual));
    std::vector<Plan> out;
    out.reserve(actual);
    for (size_t index = 0; index < actual; ++index) {
      tr_plan_t* child = nullptr;
      throw_if_error(tr_plan_partition_contiguous(handle_, count, index, &child));
      out.emplace_back(child);
    }
    return out;
  }

  Stream open_stream() const {
    tr_stream_t* stream = nullptr;
    throw_if_error(tr_stream_open_default(handle_, &stream));
    return Stream(stream);
  }

  Stream open_stream(tr_stream_config_t config) const {
    tr_stream_t* stream = nullptr;
    throw_if_error(tr_stream_open(handle_, config, &stream));
    return Stream(stream);
  }

 private:
  void reset() {
    if (handle_ != nullptr) {
      tr_plan_destroy(handle_);
      handle_ = nullptr;
    }
  }

  mutable tr_plan_t* handle_ = nullptr;
};

class Run {
 public:
  explicit Run(const std::string& path) {
    tr_run_t* handle = nullptr;
    throw_if_error(tr_run_open(path.c_str(), &handle));
    handle_ = handle;
  }

  Run(Run&& other) noexcept : handle_(std::exchange(other.handle_, nullptr)) {}
  Run& operator=(Run&& other) noexcept {
    if (this != &other) {
      reset();
      handle_ = std::exchange(other.handle_, nullptr);
    }
    return *this;
  }

  Run(const Run&) = delete;
  Run& operator=(const Run&) = delete;

  ~Run() { reset(); }

  std::vector<IsolationWindowInfo> isolation_windows() const {
    size_t count = 0;
    throw_if_error(tr_run_ms2_window_count(handle_, &count));
    std::vector<IsolationWindowInfo> out;
    out.reserve(count);
    for (size_t index = 0; index < count; ++index) {
      tr_isolation_window_info_t raw{};
      throw_if_error(tr_run_get_ms2_window(handle_, index, &raw));
      out.push_back(IsolationWindowInfo{
          raw.isolation_window_id,
          raw.window_group_id,
          raw.scan_start_index,
          raw.scan_end_index_exclusive,
          raw.isolation_target_mz,
          raw.isolation_lower_offset,
          raw.isolation_upper_offset,
      });
    }
    return out;
  }

  Plan make_plan(Representation representation, Selector selector) const {
    tr_plan_t* plan = nullptr;
    throw_if_error(tr_plan_create(
        handle_,
        static_cast<tr_representation_t>(representation),
        selector.raw,
        &plan));
    return Plan(plan);
  }

 private:
  void reset() {
    if (handle_ != nullptr) {
      tr_run_destroy(handle_);
      handle_ = nullptr;
    }
  }

  tr_run_t* handle_ = nullptr;
};

}  // namespace timsreader

#endif
