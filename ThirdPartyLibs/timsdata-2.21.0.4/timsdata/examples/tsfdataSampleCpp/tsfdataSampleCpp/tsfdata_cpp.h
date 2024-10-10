#ifndef DE_BDAL_CPP_IO_TSFDATA_CPP_H
#define DE_BDAL_CPP_IO_TSFDATA_CPP_H

/** \file
 *
 * Sample for a light-weight header-only C++ layer wrapping the C API for Bruker's TSF
 * reader DLL. You can modify this file as desired.
 *
 * See 'tsfdata.h' for more details about the underlying C API.
 *
 */

#include <stdexcept>
#include <string>
#include <cstdint>
#include <vector>
#include <limits>
#include <memory>

#include "tsfdata.h" // fundamental C API

namespace tsfdata
{
    /// \throws std::runtime_error containing the last timsdata.dll error string.
    inline void throwLastError ()
    {
        const uint32_t len = tsf_get_last_error_string(nullptr, 0);

        std::unique_ptr<char[]> buf(new char[len]);
        tsf_get_last_error_string(buf.get(), len);

        throw std::runtime_error(buf.get());
    }

    /// Reader for TSF binary data (.tsf_bin). (The SQLite file (.tsf) containing all the
    /// metadata may be opened separately using any desired SQLite API.)
    class TsfData
    {
    public:
        /// Open specified TSF analysis.
        ///
        /// \param[in] analysis_directory_name in UTF-8 encoding.
        /// \param[in] use_recalibration if DA re-calibration shall be used if available
        ///
        /// \throws std::exception in case of an error
        explicit TsfData(const std::string &analysis_directory_name, bool use_recalibration=false)
        {
            handle = tsf_open(analysis_directory_name.c_str(), use_recalibration);
            if(handle == 0)
                throwLastError();
        }

        TsfData(const TsfData&) = delete;
        TsfData operator= (const TsfData&) = delete;

        TsfData(TsfData&& other) = default;
        TsfData& operator=(TsfData&&) = default;

        /// Close TSF analysis.
        ~TsfData()
        {
            tsf_close(handle);
        }

        /// Get the C-API handle corresponding to this instance. (Caller does not get
        /// ownership of the handle.) (This call is here for the case that the user wants
        /// to call C-library functions directly.)
        uint64_t getHandle () const
        {
            return handle;
        }

        void readLineSpectrum(
            int64_t spectrum_id,
            std::vector<double>& indices,
            std::vector<float>& intensities)
        {

            // buffer-growing loop
            for (;;) {
                indices.resize(line_buffer_size);
                intensities.resize(line_buffer_size);

                const int32_t required_len = tsf_read_line_spectrum_v2(handle, spectrum_id, indices.data(), intensities.data(), static_cast<int32_t>(indices.size()));

                if (required_len < 0) {
                    throwLastError();
                }

                if (line_buffer_size >= required_len) {
                    indices.resize(required_len);
                    intensities.resize(required_len);
                    return;
                }

                line_buffer_size = required_len; // grow buffer
            }
        }

        void readLineSpectrumWithWidth(
            int64_t spectrum_id,
            std::vector<double>& indices,
            std::vector<float>& intensities,
            std::vector<float>& width)
        {
            // buffer-growing loop
            for (;;) {
                indices.resize(line_buffer_size);
                intensities.resize(line_buffer_size);
                width.resize(line_buffer_size);

                const int32_t required_len = tsf_read_line_spectrum_with_width_v2(
                    handle, spectrum_id, indices.data(), intensities.data(), width.data(), static_cast<int32_t>(indices.size()));

                if (required_len < 0) {
                    throwLastError();
                }

                if (line_buffer_size >= required_len) {
                    indices.resize(required_len);
                    intensities.resize(required_len);
                    width.resize(required_len);
                    return;
                }

                line_buffer_size = required_len; // grow buffer
            }
        }

        void readProfileSpectrum(
            int64_t spectrum_id,
            std::vector<uint32_t>& profile)
        {
            // buffer-growing loop
            for (;;) {
                profile.resize(profile_size);

                const int32_t required_len = tsf_read_profile_spectrum_v2(
                    handle, spectrum_id, profile.data(), static_cast<int32_t>(profile.size()));

                if (required_len < 0) {
                    throwLastError();
                }

                if (profile_size >= required_len) {
                    profile.resize(required_len);
                    return;
                }

                profile_size = required_len; // grow buffer
            }
        }

        #define BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(CPPNAME, CNAME) \
        void CPPNAME ( \
            int64_t frame_id,               /**< frame index */ \
            const std::vector<double> & in, /**< vector of input values (can be empty) */ \
            std::vector<double> & out )     /**< vector of corresponding output values (will be resized automatically) */ \
        { \
            doTransformation(frame_id, in, out, CNAME); \
        }

        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(indexToMz, tsf_index_to_mz)
        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(mzToIndex, tsf_mz_to_index)

    private:
        uint64_t handle {0};
        int32_t line_buffer_size{128}; // number of entries
        int32_t profile_size{0};

        void doTransformation (
            int64_t frame_id,
            const std::vector<double> & in,
            std::vector<double> & out,
            BdalTimsConversionFunction * func )
        {
            if(in.empty())
            {
                out.clear();
                return;
            }
            if(in.size() > std::numeric_limits<uint32_t>::max())
                throw std::runtime_error("Input range too large.");
            out.resize(in.size());
            func(handle, frame_id, &in[0], &out[0], uint32_t(in.size()));
        }
    };

} // namespace tsfdata

#endif // DE_BDAL_CPP_IO_TSFDATA_CPP_H
