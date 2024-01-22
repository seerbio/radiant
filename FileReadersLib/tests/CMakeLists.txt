#https://github.com/KDAB/Charm/blob/master/Tests/CMakeLists.txt for example

SET(TEST_MODULE FileReadersLib)

SET(TEST_LIBRARIES
        ${TEST_MODULE}
        Qt5::Core
        Qt5::Test
)

INCLUDE_DIRECTORIES(
        ${CMAKE_BINARY_DIR}/bin
)

#Add classes to test here.
SET(TEST_CLASSES
        CommandLineParserUtilsTests
        CSVReaderTests
        FastaReaderTests
        FeatureFinderHillReaderTests
        FragLibReaderTests
        MsCalibrationReaderTests
        MsReaderBaseTests
        MsReaderMZMLTests
        MsReaderParquetTests
        MsReaderPointerAccTests
        ParquetReaderTests
        PythiaParameterReaderTests
)

foreach(class ${TEST_CLASSES})

    SET( CLASS_SRCS ${class}.cpp)
    ADD_EXECUTABLE( ${class} ${CLASS_SRCS} )
    TARGET_LINK_LIBRARIES( ${class} ${TEST_LIBRARIES} )
    SET_TARGET_PROPERTIES(${class} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
            ${CMAKE_BINARY_DIR}/bin )
    ADD_TEST( NAME ${class} COMMAND ${class} )

endforeach()


set(FASTA_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/human_plasma_entrapment_super_trunc.fasta")
configure_file(${FASTA_PATH} ${CMAKE_BINARY_DIR}/bin/human_plasma_entrapment_super_trunc.fasta COPYONLY)

set(PYTHIA_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/test_params_wide_window.pythiaConfig")
configure_file(${PYTHIA_PATH} ${CMAKE_BINARY_DIR}/bin/test_params_wide_window.pythiaConfig COPYONLY)

set(MS_PARQUET "${PROJECT_SOURCE_DIR}/tests/TestFiles/EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF")
configure_file(${MS_PARQUET} ${CMAKE_BINARY_DIR}/bin/EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF COPYONLY)

set(FRAG_LIB_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/FragLibReaderTests.fragLibFF")
configure_file(${FRAG_LIB_PATH} ${CMAKE_BINARY_DIR}/bin/FragLibReaderTests.fragLibFF COPYONLY)

set(MZML_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/1min.mzML")
configure_file(${MZML_PATH} ${CMAKE_BINARY_DIR}/bin/1min.mzML COPYONLY)

set(FRAME_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/framePoints.frame")
configure_file(${FRAME_PATH} ${CMAKE_BINARY_DIR}/bin/framePoints.frame COPYONLY)
