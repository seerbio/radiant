#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "UtilsLib" for configuration "RelWithDebInfo"
set_property(TARGET UtilsLib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(UtilsLib PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libUtilsLib.so"
  IMPORTED_SONAME_RELWITHDEBINFO "libUtilsLib.so"
  )

list(APPEND _cmake_import_check_targets UtilsLib )
list(APPEND _cmake_import_check_files_for_UtilsLib "${_IMPORT_PREFIX}/lib/libUtilsLib.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
