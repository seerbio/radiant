#https://github.com/KDAB/Charm/blob/master/Tests/CMakeLists.txt for example

SET(TEST_MODULE AlgorithmsLib)

SET(TEST_LIBRARIES
        ${TEST_MODULE}
        Qt5::Core
        Qt5::Test
)

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/bin)

#Add classes to test here.
SET(TEST_CLASSES
        BiophysicalCalcsTests
        FeatureFinderHillBuilderTests
        FragLibraryTronTests
        FragmentLibraryRTreeTests
        MsScansDenoiseTronTests
        PeptidesLibraryTronTests
        ProteinDigestomaticTests
)

foreach(class ${TEST_CLASSES})

    SET( CLASS_SRCS ${class}.cpp)
    ADD_EXECUTABLE( ${class} ${CLASS_SRCS} )
    TARGET_LINK_LIBRARIES( ${class} ${TEST_LIBRARIES} )
    SET_TARGET_PROPERTIES(${class} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
            ${CMAKE_BINARY_DIR}/bin )
    ADD_TEST( NAME ${class} COMMAND ${class} )

endforeach()
