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
        CSVReaderTests
        FastaReaderTests
        FeatureFinderHillReaderTests
        FragLibReaderTests
        JsonParameterReaderTests
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

set(FASTA_PEPTIDES_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/human_plasma_entrapment_super_trunc.fasta.csv")
configure_file(${FASTA_PEPTIDES_PATH} ${CMAKE_BINARY_DIR}/bin/human_plasma_entrapment_super_trunc.fasta.csv COPYONLY)

set(FASTA_PEPLIB_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/human_plasma_entrapment_super_trunc.fasta.pepLib")
configure_file(${FASTA_PEPLIB_PATH} ${CMAKE_BINARY_DIR}/bin/human_plasma_entrapment_super_trunc.fasta.pepLib COPYONLY)

set(FASTA_FRAGLIB_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/human_plasma_entrapment_super_trunc.fasta.csv.fragLib")
configure_file(${FASTA_FRAGLIB_PATH} ${CMAKE_BINARY_DIR}/bin/human_plasma_entrapment_super_trunc.fasta.csv.fragLib COPYONLY)

set(JSON_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/DigestParamsTest.dgst")
configure_file(${JSON_PATH} ${CMAKE_BINARY_DIR}/bin/DigestParamsTest.dgst COPYONLY)

set(PYTHIA_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/DigestParamsTest2.pythia")
configure_file(${PYTHIA_PATH} ${CMAKE_BINARY_DIR}/bin/DigestParamsTest2.pythia COPYONLY)

set(MS_PARQUET "${PROJECT_SOURCE_DIR}/tests/TestFiles/EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF")
configure_file(${MS_PARQUET} ${CMAKE_BINARY_DIR}/bin/EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF COPYONLY)

set(FRAG_LIB_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/FragLibReaderTests.fragLibFF")
configure_file(${FRAG_LIB_PATH} ${CMAKE_BINARY_DIR}/bin/FragLibReaderTests.fragLibFF COPYONLY)

set(TEST_PARAMS_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/WorkFlowTestsParams.pythia")
configure_file(${TEST_PARAMS_PATH} ${CMAKE_BINARY_DIR}/bin/WorkFlowTestsParams.pythia COPYONLY)

set(PEPTIDES_CSV_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/peptides.csv")
configure_file(${PEPTIDES_CSV_PATH} ${CMAKE_BINARY_DIR}/bin/peptides.csv COPYONLY)

set(TARGETS_MZ_CE_CSV_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/target_mz_ce.csv")
configure_file(${TARGETS_MZ_CE_CSV_PATH} ${CMAKE_BINARY_DIR}/bin/target_mz_ce.csv COPYONLY)

set(FRAME_SCANS_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScans")
configure_file(${FRAME_SCANS_PATH} ${CMAKE_BINARY_DIR}/bin/EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScans COPYONLY)

set(FRAME_SCORES_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScores")
configure_file(${FRAME_SCORES_PATH} ${CMAKE_BINARY_DIR}/bin/EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScores COPYONLY)

set(PEP_LIB_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/2022_02_22_Homo_sapiens_UP000005640.fasta.pepLib")
configure_file(${PEP_LIB_PATH} ${CMAKE_BINARY_DIR}/bin/2022_02_22_Homo_sapiens_UP000005640.fasta.pepLib COPYONLY)

