cmake_minimum_required(VERSION 3.18)
PROJECT(MachineLrnLib
        VERSION "1.0"
        )

INCLUDE(GenerateExportHeader)
INCLUDE(CMakePackageConfigHelpers)

SET(CMAKE_CXX_STANDARD 14)
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
        src/NearestNeighborsSearch.cpp
        src/NeuralNetModel.cpp
        )

SET(HEADER_FILES
        src/NearestNeighborsSearch.h
        src/NeuralNetModel.h
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
        "${EIGEN_DIRECTORY}"
        "${FUNCTIONAL_PLUS_DIRECTORY}"
        "${FRUGALLY_DEEP_DIRECTORY}"
        "${JSON_DIRECTORY}"
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        "${NANOFLANN_DIRECTORY}"
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
        EigenLib
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

# for MSVS
SOURCE_GROUP("Source Files" FILES ${${PROJECT_NAME}_LIB_SOURCES})

CONFIGURE_FILE(${MODELS_DIRECTORY}/rnn_linear_charge_w_precursors_nce_1.hdf5.json ${CMAKE_BINARY_DIR}/bin/rnn_linear_charge_w_precursors_nce_1.hdf5.json  COPYONLY)
CONFIGURE_FILE(${MODELS_DIRECTORY}/rnn_linear_charge_w_precursors_nce_2.hdf5.json ${CMAKE_BINARY_DIR}/bin/rnn_linear_charge_w_precursors_nce_2.hdf5.json  COPYONLY)
CONFIGURE_FILE(${MODELS_DIRECTORY}/rnn_linear_charge_w_precursors_nce_3.hdf5.json ${CMAKE_BINARY_DIR}/bin/rnn_linear_charge_w_precursors_nce_3.hdf5.json  COPYONLY)
CONFIGURE_FILE(${MODELS_DIRECTORY}/rnn_linear_charge_w_precursors_nce_4.hdf5.json ${CMAKE_BINARY_DIR}/bin/rnn_linear_charge_w_precursors_nce_4.hdf5.json  COPYONLY)
CONFIGURE_FILE(${MODELS_DIRECTORY}/MS2_Mono_Model_2.json ${CMAKE_BINARY_DIR}/bin/MS2_Mono_Model.json  COPYONLY)
CONFIGURE_FILE(${MODELS_DIRECTORY}/MS2_Charge_Model_4.json ${CMAKE_BINARY_DIR}/bin/MS2_Charge_Model.json  COPYONLY)
CONFIGURE_FILE(${MODELS_DIRECTORY}/iRT_Model.json ${CMAKE_BINARY_DIR}/bin/iRT_Model.json  COPYONLY)

ADD_SUBDIRECTORY("tests")
