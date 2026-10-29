ADD_EXECUTABLE(RadiantCommandLineParserTests
        CommandLineParserTests.cpp
        ../src/CommandLineParser.cpp
        )

TARGET_INCLUDE_DIRECTORIES(RadiantCommandLineParserTests PRIVATE
        ../src
        )

TARGET_LINK_LIBRARIES(RadiantCommandLineParserTests
        Qt5::Core
        Qt5::Test
        FileReadersLib
        UtilsLib
        )

SET_TARGET_PROPERTIES(RadiantCommandLineParserTests PROPERTIES
        AUTOMOC ON
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )

ADD_TEST(NAME RadiantCommandLineParserTests COMMAND RadiantCommandLineParserTests)

FIND_PACKAGE(Python3 COMPONENTS Interpreter QUIET)
IF(Python3_Interpreter_FOUND)
    ADD_TEST(NAME RadiantBatchCliTests
            COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/BatchCliTests.py
            --binary $<TARGET_FILE:RadiantDIA>
            --fixtures ${CMAKE_BINARY_DIR}/bin
            )
ENDIF()
