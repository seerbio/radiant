# -*- coding: utf-8 -*-
"""Python wrapper for timsdata.dll for reading tsf"""

import numpy as np
import sqlite3
import os, sys
from ctypes import *
from pathlib import Path

if sys.platform[:5] == "win32":
    libname = "timsdata.dll"
elif sys.platform[:5] == "linux":
    libname = "libtimsdata.so"
else:
    raise Exception("Unsupported platform.")

path = Path(Path(__file__).parent.absolute(), 'libs', libname).as_posix()
if os.path.exists(path):
    dll = cdll.LoadLibrary(path)
else:
    dll = cdll.LoadLibrary(libname)

dll.tsf_open.argtypes = [ c_char_p, c_uint32 ]
dll.tsf_open.restype = c_uint64
dll.tsf_close.argtypes = [ c_uint64 ]
dll.tsf_close.restype = None
dll.tsf_get_last_error_string.argtypes = [ c_char_p, c_uint32 ]
dll.tsf_get_last_error_string.restype = c_uint32
dll.tsf_has_recalibrated_state.argtypes = [ c_uint64 ]
dll.tsf_has_recalibrated_state.restype = c_uint32
dll.tsf_read_line_spectrum_v2.argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_float), c_int32 ]
dll.tsf_read_line_spectrum_v2.restype = c_int32
dll.tsf_read_line_spectrum_with_width_v2.argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_float), POINTER(c_float), c_int32 ]
dll.tsf_read_line_spectrum_with_width_v2.restype = c_int32
dll.tsf_read_profile_spectrum_v2.argtypes = [ c_uint64, c_int64, POINTER(c_uint32), c_int32 ]
dll.tsf_read_profile_spectrum_v2.restype = c_int32

convfunc_argtypes = [ c_uint64, c_int64, POINTER(c_double), POINTER(c_double), c_uint32 ]

dll.tsf_index_to_mz.argtypes = convfunc_argtypes
dll.tsf_index_to_mz.restype = c_uint32
dll.tsf_mz_to_index.argtypes = convfunc_argtypes
dll.tsf_mz_to_index.restype = c_uint32

def _throwLastTsfDataError (dll_handle):
    """Throw last TimsData error string as an exception."""

    len = dll_handle.tsf_get_last_error_string(None, 0)
    buf = create_string_buffer(len)
    dll_handle.tsf_get_last_error_string(buf, len)
    raise RuntimeError(buf.value)

class TsfData:

    def __init__ (self, analysis_directory, use_recalibrated_state=False):

        if sys.version_info.major == 2:
            if not isinstance(analysis_directory, unicode):
                raise ValueError("analysis_directory must be a Unicode string.")
        if sys.version_info.major == 3:
            if not isinstance(analysis_directory, str):
                raise ValueError("analysis_directory must be a string.")

        self.dll = dll

        self.handle = self.dll.tsf_open(
            analysis_directory.encode('utf-8'),
            1 if use_recalibrated_state else 0 )
        if self.handle == 0:
            _throwLastTsfDataError(self.dll)

        self.conn = sqlite3.connect(os.path.join(analysis_directory, "analysis.tsf"))

        self.line_buffer_size = 1024 # may grow in read...Spectrum()
        self.profile_buffer_size = 1024 # may grow in read...Spectrum()

    def __enter__(self):
        return self
        
    def __exit__(self, type, value, traceback):
        self.close()
        
    def __del__ (self):
        self.close()
        
    def close(self):
        if hasattr(self, 'handle') and self.handle is not None:
            self.dll.tsf_close(self.handle)
            self.handle = None
        if hasattr(self, 'conn') and self.conn is not None:
            self.conn.close()
            self.conn = None

    def __callConversionFunc (self, frame_id, input_data, func):

        if type(input_data) is np.ndarray and input_data.dtype == np.float64:
            # already "native" format understood by DLL -> avoid extra copy
            in_array = input_data
        else:
            # convert data to format understood by DLL:
            in_array = np.array(input_data, dtype=np.float64)

        cnt = len(in_array)
        out = np.empty(shape=cnt, dtype=np.float64)
        success = func(self.handle, frame_id,
                       in_array.ctypes.data_as(POINTER(c_double)),
                       out.ctypes.data_as(POINTER(c_double)),
                       cnt)

        if success == 0:
            _throwLastTsfDataError(self.dll)

        return out

    def indexToMz (self, frame_id, indices):
        return self.__callConversionFunc(frame_id, indices, self.dll.tsf_index_to_mz)

    def mzToIndex (self, frame_id, mzs):
        return self.__callConversionFunc(frame_id, mzs, self.dll.tsf_mz_to_index)

    # Output: tuple of lists (indices, intensities)
    def readLineSpectrum (self, frame_id):
        # buffer-growing loop
        while True:
            cnt = int(self.profile_buffer_size) # necessary cast to run with python 3.5
            index_buf = np.empty(shape=cnt, dtype=np.float64)
            intensity_buf = np.empty(shape=cnt, dtype=np.float32)

            required_len = self.dll.tsf_read_line_spectrum_v2(self.handle, frame_id, index_buf.ctypes.data_as(POINTER(c_double)), intensity_buf.ctypes.data_as(POINTER(c_float)), self.profile_buffer_size)

            if required_len < 0:
                _throwLastTsfDataError(self.dll)

            if required_len > self.profile_buffer_size:
                if required_len > 16777216:
                    # arbitrary limit for now...
                    raise RuntimeError("Maximum expected frame size exceeded.")
                self.profile_buffer_size = required_len # grow buffer
            else:
                break

        return (index_buf[0 : required_len], intensity_buf[0 : required_len])

        # Output: tuple of lists (indices, intensities, widths)
    def readLineSpectrumWithWidth (self, frame_id):
        # buffer-growing loop
        while True:
            cnt = int(self.profile_buffer_size) # necessary cast to run with python 3.5
            index_buf = np.empty(shape=cnt, dtype=np.float64)
            intensity_buf = np.empty(shape=cnt, dtype=np.float32)
            width_buf = np.empty(shape=cnt, dtype=np.float32)

            required_len = self.dll.tsf_read_line_spectrum_with_width_v2(
                self.handle, frame_id, index_buf.ctypes.data_as(POINTER(c_double)), intensity_buf.ctypes.data_as(POINTER(c_float)), width_buf.ctypes.data_as(POINTER(c_float)), self.profile_buffer_size)

            if required_len < 0:
                _throwLastTsfDataError(self.dll)

            if required_len > self.profile_buffer_size:
                if required_len > 16777216:
                    # arbitrary limit for now...
                    raise RuntimeError("Maximum expected frame size exceeded.")
                self.profile_buffer_size = required_len # grow buffer
            else:
                break

        return (index_buf[0 : required_len], intensity_buf[0 : required_len], width_buf[0 : required_len])

    # Output: tuple of lists (indices, intensities)
    def readProfileSpectrum (self, frame_id):
        # buffer-growing loop
        while True:
            cnt = int(self.profile_buffer_size) # necessary cast to run with python 3.5
            intensity_buf = np.empty(shape=cnt, dtype=np.uint32)

            required_len = self.dll.tsf_read_profile_spectrum_v2(self.handle, frame_id, intensity_buf.ctypes.data_as(POINTER(c_uint32)), self.profile_buffer_size)

            if required_len < 0:
                _throwLastTsfDataError(self.dll)

            if required_len > self.profile_buffer_size:
                if required_len > 16777216:
                    # arbitrary limit for now...
                    raise RuntimeError("Maximum expected frame size exceeded.")
                self.profile_buffer_size = required_len # grow buffer
            else:
                break

        return intensity_buf[0 : required_len]
