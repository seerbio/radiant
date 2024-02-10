cmake_minimum_required(VERSION 3.18)
PROJECT(AlgorithmsFFLib
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
        src/BiophysicalCalcs.cpp
        src/CandidateScores.cpp
        src/CandidateScorertron.cpp
        src/ClassifierWeightsManager.cpp
        src/DiscriminantScoretron.cpp
        src/DeMultiplexScanerator.cpp
        src/DeMultiplexScaneratorManager.cpp
        src/FDRCLassifierNeuralNet.cpp
        src/FeatureFinderHill.cpp
        src/FeatureFinderHillBuilder.cpp
        src/InterferingCandidatesEliminatomatic.cpp
        src/IsotopicDistributionBuilder.cpp
        src/MS2Ion.cpp
        src/MsCalibratomatic.cpp
        src/MsFrame.cpp
        src/MsUtils.cpp
        src/PeakIntegratomatic.cpp
        src/ProteinDigestomatic.cpp
        src/QValueSettertron.cpp
        src/ScoreOverseer.cpp
        src/TargetDecoyCandidatePair.cpp
        src/TargetDecoyCandidatePairManager.cpp
        src/TargetDecoyCandidatePairScoretron.cpp
        src/TurboXIC.cpp
        src/XYMappermatic.cpp
)

SET(HEADER_FILES
        src/BiophysicalCalcs.h
        src/CandidateScores.h
        src/CandidateScorertron.h
        src/ClassifierWeightsManager.h
        src/DeMultiplexScanerator.h
        src/DeMultiplexScaneratorManager.h
        src/DiscriminantScoretron.h
        src/FDRCLassifierNeuralNet.h
        src/FeatureFinderHill.h
        src/FeatureFinderHillBuilder.h
        src/InterferingCandidatesEliminatomatic.h
        src/IsotopicDistributionBuilder.h
        src/MS2Ion.h
        src/MsCalibratomatic.h
        src/MsFrame.h
        src/MsUtils.h
        src/PeakIntegratomatic.h
        src/ProteinDigestomatic.h
        src/QValueSettertron.h
        src/ScoreOverseer.cpp
        src/TargetDecoyCandidatePair.h
        src/TargetDecoyCandidatePairManager.h
        src/TargetDecoyCandidatePairScoretron.h
        src/TurboXIC.h
        src/XYMappermatic.h
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
        )

FIND_PACKAGE(Qt5 COMPONENTS
        Concurrent
        Core
        Test
        REQUIRED
)


ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_LIB_SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
        Qt5::Concurrent
        Qt5::Core
        Qt5::Test
        ChemLib
        EigenLib
        FileReadersLib
        PyTorchLib
        UtilsLib
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
