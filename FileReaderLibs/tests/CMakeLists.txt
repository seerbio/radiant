#https://github.com/KDAB/Charm/blob/master/Tests/CMakeLists.txt for example

SET(TEST_MODULE FileReadersLib)

SET(TEST_LIBRARIES
        ${TEST_MODULE}
        Qt5::Core
        Qt5::Test
)

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/bin)

#Add classes to test here.
SET(TEST_CLASSES
        FastaReaderTests
        JsonParameterReaderTests
        MsParquetReaderTests
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

set(JSON_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/DigestParamsTest.dgst")
configure_file(${JSON_PATH} ${CMAKE_BINARY_DIR}/bin/DigestParamsTest.dgst COPYONLY)

set(PYTHIA_PATH "${PROJECT_SOURCE_DIR}/tests/TestFiles/DigestParamsTest2.pythia")
configure_file(${PYTHIA_PATH} ${CMAKE_BINARY_DIR}/bin/DigestParamsTest2.pythia COPYONLY)
