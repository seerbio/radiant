cmake_minimum_required(VERSION 3.18)
PROJECT(FileReadersLib
        VERSION "1.0"
        )

INCLUDE(GenerateExportHeader)
INCLUDE(CMakePackageConfigHelpers)

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTORCC ON)
SET(CMAKE_AUTOUIC ON)

SET(CMAKE_AUTOMOC_MACRO_NAMES
        "Q_OBJECT"
        "Q_GADGET"
        "Q_NAMESPACE"
        "Q_NAMESPACE_EXPORT"
        )

SET(SOURCE_FILES

        src/SpecLibSrc/Library.cpp
        src/SpecLibSrc/Readers.cpp
        src/SpecLibSrc/SpecLibStructures.cpp

        src/MsCalibrationReader.cpp
        src/CommandLineParserUtils.cpp
        src/CSVReader.cpp
        src/FastaReader.cpp
        src/FeatureFinderHillReader.cpp
        src/FragLibReader.cpp
        src/FragLibReaderRow.cpp
        src/FragLibTsvReader.cpp
        src/MsReaderParquet.cpp
        src/MsReaderBase.cpp
#        src/MsReaderBrukerTims.cpp
        src/MsReaderMzML.cpp
        src/MsReaderMzMLLazyLoad.cpp
        src/MsReaderMzMLMapped.cpp
        src/MsReaderPointerAcc.cpp
        src/ParquetReader.cpp
        src/PythiaParameterReader.cpp
        src/QuanReader.cpp
        src/SpecLibReader.cpp
        ${TIMS_CPPSQL_SOURCE}
        )

SET(HEADER_FILES

        src/SpecLibSrc/Library.h
        src/SpecLibSrc/Readers.h
        src/SpecLibSrc/SpecLibStructures.h

        src/MsCalibrationReader.h
        src/CommandLineParserUtils.h
        src/CSVReader.h
        src/FastaReader.h
        src/FeatureFinderHillReader.h
        src/FragLibReader.h
        src/FragLibReaderRow.h
        src/FragLibTsvReader.h
        src/MsReaderParquet.h
        src/MsReaderBase.h
#        src/MsReaderBrukerTims.h
        src/MsReaderMzML.h
        src/MsReaderMzMLLazyLoad.h
        src/MsReaderMzMLMapped.h
        src/MsReaderPointerAcc.h
        src/ParquetReader.h
        src/PythiaParameterReader.h
        src/QuanReader.h
        src/SpecLibReader.h
        src/toml.hpp
#        ${TIMS_CPPSQL_HEADER}
        )

SET(${PROJECT_NAME}_LIB_SOURCES
        ${SOURCE_FILES}
        ${HEADER_FILES}
)

SET(DLL_LIB_INCLUDE_DIRS
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

INCLUDE_DIRECTORIES(
        DLL_LIB_INCLUDE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        ${EIGEN_DIRECTORY}
        ${NANOFLANN_DIRECTORY}
#        ${TIMS_INCLUDE_DIRECTORY}
#        ${TIMS_CPPSQL_DIRECTORY}
        ${LIBXML2_INCLUDE_DIR}
)

FIND_PACKAGE(Parquet REQUIRED)
FIND_PACKAGE(LibXml2 REQUIRED)
FIND_PACKAGE(Qt5 COMPONENTS
        Concurrent
        Test
        REQUIRED
)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_LIB_SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
        Qt5::Concurrent
        ChemLib
        EigenLib
        UtilsLib
        Parquet::parquet_static
#        ${TIMS_LIB}
        ${LIBXML2_LIBRARIES}
        sqlite3
        )

# This puts binaries and exes in bin directory.  Otherwise won't run
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/bin
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )

SET_PROPERTY(TARGET ${PROJECT_NAME} PROPERTY FOLDER "Libs")

STRING(TOUPPER ${PROJECT_NAME}_EXPORTS EXPORT_NAME)
GENERATE_EXPORT_HEADER(${PROJECT_NAME}
        BASE_NAME ${PROJECT_NAME}
        EXPORT_MACRO_NAME ${EXPORT_NAME}
        EXPORT_FILE_NAME ${PROJECT_NAME}_Exports.h
        STATIC_DEFINE ${PROJECT_NAME}_EXPORTS_BUILT_AS_STATIC)

TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME}
        PUBLIC "$<BUILD_INTERFACE:${DLL_LIB_INCLUDE_DIRS}/src>"
        INTERFACE "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${ROOT_PROJECT_NAME}>"
        )

SET(CMAKECONFIG_INSTALL_PREFIX "lib/cmake/${PROJECT_NAME}")
SET(CMAKECONFIG_INSTALL_DIR "${CMAKECONFIG_INSTALL_PREFIX}")


CONFIGURE_PACKAGE_CONFIG_FILE(${PROJECT_NAME}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR}
        )

WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY ExactVersion
        )


INSTALL(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR}
        COMPONENT Devel
        )

INSTALL(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        ${INSTALL_TARGETS_DEFAULT_ARGS}
        )

INSTALL(EXPORT ${PROJECT_NAME}Targets
        DESTINATION "${CMAKECONFIG_INSTALL_DIR}"
        FILE ${PROJECT_NAME}Targets.cmake
        COMPONENT Devel
        )

ADD_SUBDIRECTORY("tests")
