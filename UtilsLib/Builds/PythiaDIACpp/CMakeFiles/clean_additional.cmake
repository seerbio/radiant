# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "RelWithDebInfo")
  file(REMOVE_RECURSE
  "CMakeFiles/UtilsLib_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/UtilsLib_autogen.dir/ParseCache.txt"
  "UtilsLib_autogen"
  "tests/CMakeFiles/ErrorUtilsTests_autogen.dir/AutogenUsed.txt"
  "tests/CMakeFiles/ErrorUtilsTests_autogen.dir/ParseCache.txt"
  "tests/CMakeFiles/MathUtilsTests_autogen.dir/AutogenUsed.txt"
  "tests/CMakeFiles/MathUtilsTests_autogen.dir/ParseCache.txt"
  "tests/CMakeFiles/SqlUtilsTests_autogen.dir/AutogenUsed.txt"
  "tests/CMakeFiles/SqlUtilsTests_autogen.dir/ParseCache.txt"
  "tests/ErrorUtilsTests_autogen"
  "tests/MathUtilsTests_autogen"
  "tests/SqlUtilsTests_autogen"
  )
endif()
