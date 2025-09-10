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
#        src/CandidateScores.cpp
        src/CandidateScoresDDA.cpp
        src/CandidateScoresV2.cpp
        src/CandidateScoresFeatureManager.cpp
#        src/CandidateScorertron.cpp
        src/CandidateScorertronV2.cpp
#        src/Centroidotron.cpp
#        src/ClassifierWeightsManager.cpp
#        src/Deconvolvotron.cpp
#        src/DeIsotopotron.cpp
#        src/DiscriminantScoretron.cpp
#        src/DeMultiplexScanerator.cpp
#        src/DeMultiplexScaneratorManager.cpp
#        src/FDRCLassifierNeuralNet.cpp
        src/FeatureFinderHill.cpp
        src/FeatureFinderHillBuilder.cpp
#        src/IonMobilitron.cpp
#        src/InterferingCandidatesEliminatomatic.cpp
#        src/IsotopicDistributionBuilder.cpp
        src/MS2Ion.cpp
        src/Ms2IonFraggertronManager.cpp
        src/MsCalibratomatic.cpp
        src/MsFraggertron.cpp
#        src/MsFrame.cpp
        src/MsFrameV2.cpp
        src/MsUtils.cpp
        src/PeakIntegratomatic.cpp
        src/PrecursorMzArrangetron.cpp
#        src/QuanFileBuilder.cpp
#        src/QuanTransitionRefinertron.cpp
#        src/QValueSettertron.cpp
        src/QValueSettertronV2.cpp
#        src/SequenceSubstringSearchomatic.cpp
#        src/SpectrumCentricMzTargetFrameSearch.cpp
        src/TargetDecoyCandidatePair.cpp
        src/TargetDecoyCandidatePairManager.cpp
#        src/TargetDecoyCandidatePairScoretron.cpp
        src/TargetDecoyCandidatePairScoretronV2.cpp
        src/TurboXIC.cpp
#        src/TurboXIC2D.cpp
#        src/XICPeakManager.cpp
        src/XYMappermatic.cpp
)

SET(HEADER_FILES
#        src/CandidateScores.h
        src/CandidateScoresDDA.h
        src/CandidateScoresV2.h
        src/CandidateScoresFeatureManager.h
#        src/CandidateScorertron.h
        src/CandidateScorertronV2.h
#        src/Centroidotron.h
#        src/ClassifierWeightsManager.h
#        src/Deconvolvotron.h
#        src/DeIsotopotron.h
#        src/DeMultiplexScanerator.h
#        src/DeMultiplexScaneratorManager.h
#        src/DiscriminantScoretron.h
#        src/FDRCLassifierNeuralNet.h
        src/FeatureFinderHill.h
        src/FeatureFinderHillBuilder.h
#        src/InterferingCandidatesEliminatomatic.h
#        src/IonMobilitron.h
#        src/IsotopicDistributionBuilder.h
        src/MS2Ion.h
        src/Ms2IonFraggertronManager.h
        src/MsCalibratomatic.h
        src/MsFraggertron.h
#        src/MsFrame.h
        src/MsFrameV2.h
        src/MsUtils.h
        src/PeakIntegratomatic.h
        src/PrecursorMzArrangetron.h
#        src/QuanFileBuilder.h
#        src/QuanTransitionRefinertron.h
#        src/QValueSettertron.h
        src/QValueSettertronV2.h
#        src/SequenceSubstringSearchomatic.h
#        src/SpectrumCentricMzTargetFrameSearch.h
        src/TargetDecoyCandidatePair.h
        src/TargetDecoyCandidatePairManager.h
#        src/TargetDecoyCandidatePairScoretron.h
        src/TargetDecoyCandidatePairScoretronV2.h
        src/TurboXIC.h
#        src/TurboXIC2D.h
#        src/XICPeakManager.h
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

FIND_PACKAGE (SeqAn REQUIRED)

FIND_PACKAGE(Qt5 COMPONENTS
        Concurrent
        Core
        Test
        REQUIRED
)

INCLUDE_DIRECTORIES(
        DLL_LIB_INCLUDE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        ${EIGEN_DIRECTORY}
#        ${NANOFLANN_DIRECTORY}
#        ${SEQAN_INCLUDE_DIRS}
        )

ADD_DEFINITIONS (${SEQAN_DEFINITIONS})
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_LIB_SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
        Qt5::Concurrent
        Qt5::Core
        Qt5::Test
        AVXLib
        ChemLib
        EigenLib
        FileReadersLib
#        PyTorchLib
        UtilsLib
#        ${SEQAN_LIBRARIES}
        )

# This puts binaries and exes in bin directory.  Otherwise won't run
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/bin
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        )

include(CheckCXXCompilerFlag)
# Check for AVX
check_cxx_compiler_flag("-mavx" COMPILER_SUPPORTS_AVX)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)

if (COMPILER_SUPPORTS_AVX512)
    message(STATUS "Enabling AVX-512, AVX2, and AVX support")
    add_compile_definitions(AVX512=1)
    target_compile_options(${PROJECT_NAME} PRIVATE -mfma -mavx2 -mavx)
elseif (COMPILER_SUPPORTS_AVX2)
    message(STATUS "Enabling AVX2 and AVX support")
    target_compile_options(${PROJECT_NAME} PRIVATE -mavx2 -mfma -mavx)
elseif (COMPILER_SUPPORTS_AVX)
    message(STATUS "Enabling AVX support")
    target_compile_options(${PROJECT_NAME} PRIVATE -mavx)
elseif (MSVC)
    message(STATUS "Enabling AVX support for MSVC")
    target_compile_options(${PROJECT_NAME} PRIVATE /arch:AVX)
else()
    message(WARNING "AVX, AVX2, and AVX-512 are not supported by your compiler")
endif()


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
