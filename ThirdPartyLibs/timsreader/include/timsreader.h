#ifndef TIMSREADER_H
#define TIMSREADER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ArrowArray {
  int64_t length;
  int64_t null_count;
  int64_t offset;
  int64_t n_buffers;
  int64_t n_children;
  const void** buffers;
  struct ArrowArray** children;
  struct ArrowArray* dictionary;
  void (*release)(struct ArrowArray*);
  void* private_data;
};

struct ArrowSchema {
  const char* format;
  const char* name;
  const char* metadata;
  int64_t flags;
  int64_t n_children;
  struct ArrowSchema** children;
  struct ArrowSchema* dictionary;
  void (*release)(struct ArrowSchema*);
  void* private_data;
};

typedef enum tr_status_t {
  TR_STATUS_OK = 0,
  TR_STATUS_END = 1,
  TR_STATUS_INVALID_ARGUMENT = 2,
  TR_STATUS_ERROR = 3
} tr_status_t;

typedef enum tr_representation_t {
  TR_DIRECT_SUBSCANS = 0,
  TR_FLATTENED_WINDOW_SPECTRUM = 1,
  TR_CENTROIDED_FRAME_SPECTRUM = 2
} tr_representation_t;

typedef enum tr_selector_kind_t {
  TR_SELECTOR_MS1 = 0,
  /* Broad MS2 selector for processed representations. */
  TR_SELECTOR_MS2 = 1,
  /* Logical-window selector retained for direct-mode reading. */
  TR_SELECTOR_MS2_WINDOW = 2
} tr_selector_kind_t;

typedef enum tr_mz_tolerance_kind_t {
  TR_MZ_TOLERANCE_PPM = 0,
  TR_MZ_TOLERANCE_BINS = 1
} tr_mz_tolerance_kind_t;

typedef struct tr_selector_t {
  tr_selector_kind_t kind;
  /* Used only when kind == TR_SELECTOR_MS2_WINDOW. */
  uint32_t isolation_window_id;
} tr_selector_t;

typedef struct tr_stream_config_t {
  size_t batch_row_count;
  size_t centroid_max_peaks;
  tr_mz_tolerance_kind_t centroid_mz_tolerance_kind;
  double centroid_mz_tolerance_value;
  double centroid_im_pct_tol;
  uint32_t centroid_early_stop_iterations;
  bool centroid_window_cap_enabled;
  uint32_t centroid_window_cap_max_peaks;
  float centroid_window_cap_window_da;
} tr_stream_config_t;

typedef struct tr_isolation_window_info_t {
  uint32_t isolation_window_id;
  int32_t window_group_id;
  uint16_t scan_start_index;
  uint16_t scan_end_index_exclusive;
  float isolation_target_mz;
  float isolation_lower_offset;
  float isolation_upper_offset;
} tr_isolation_window_info_t;

typedef struct tr_run_t tr_run_t;
typedef struct tr_plan_t tr_plan_t;
typedef struct tr_stream_t tr_stream_t;
typedef struct tr_owned_batch_t tr_owned_batch_t;

const char* tr_last_error_message(void);
tr_stream_config_t tr_stream_config_default(void);
tr_stream_config_t tr_stream_config_default_for_selector(tr_selector_t selector);

tr_status_t tr_run_open(const char* path, tr_run_t** out_run);
void tr_run_destroy(tr_run_t* run);
tr_status_t tr_run_ms2_window_count(tr_run_t* run, size_t* out_count);
tr_status_t tr_run_get_ms2_window(
    tr_run_t* run,
    size_t index,
    tr_isolation_window_info_t* out_info);

tr_status_t tr_plan_create(
    tr_run_t* run,
    tr_representation_t representation,
    tr_selector_t selector,
    tr_plan_t** out_plan);
void tr_plan_destroy(tr_plan_t* plan);
tr_status_t tr_plan_row_count(tr_plan_t* plan, size_t* out_count);
tr_status_t tr_plan_default_stream_config(
    tr_plan_t* plan,
    tr_stream_config_t* out_config);
/* Return the set of logical MS2 window ids that may appear in this partition. */
tr_status_t tr_plan_isolation_window_count(
    tr_plan_t* plan,
    size_t* out_count);
tr_status_t tr_plan_get_isolation_window_id(
    tr_plan_t* plan,
    size_t index,
    uint32_t* out_window_id);
/* Split a plan into at most requested_count frame-aligned partitions. */
tr_status_t tr_plan_partition_count(
    tr_plan_t* plan,
    size_t requested_count,
    size_t* out_count);
tr_status_t tr_plan_partition_contiguous(
    tr_plan_t* plan,
    size_t requested_count,
    size_t partition_index,
    tr_plan_t** out_partition);

tr_status_t tr_stream_open(
    tr_plan_t* plan,
    tr_stream_config_t config,
    tr_stream_t** out_stream);
tr_status_t tr_stream_open_default(
    tr_plan_t* plan,
    tr_stream_t** out_stream);
void tr_stream_destroy(tr_stream_t* stream);
tr_status_t tr_stream_next_borrowed(
    tr_stream_t* stream,
    const struct ArrowArray** out_array,
    const struct ArrowSchema** out_schema);
tr_status_t tr_stream_next_owned(
    tr_stream_t* stream,
    tr_owned_batch_t** out_batch);

void tr_owned_batch_destroy(tr_owned_batch_t* batch);
tr_status_t tr_owned_batch_export(
    tr_owned_batch_t* batch,
    const struct ArrowArray** out_array,
    const struct ArrowSchema** out_schema);

#ifdef __cplusplus
}
#endif

#endif
