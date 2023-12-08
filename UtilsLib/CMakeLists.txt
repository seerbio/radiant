cmake_minimum_required(VERSION 3.18)
PROJECT(UtilsLib
        VERSION "1.0"
        )

INCLUDE(GenerateExportHeader)
INCLUDE(CMakePackageConfigHelpers)

enable_testing(true)

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
        src/Error.cpp
        src/ErrorUtils.cpp
        src/GlobalSettings.cpp
        src/MathUtils.cpp
        src/ParallelUtils.cpp
        src/PointDF.cpp
        src/SqlUtils.cpp
        src/StringUtils.cpp
        )

SET(HEADER_FILES
        src/Error.h
        src/ErrorUtils.h
        src/GlobalSettings.h
        src/MathUtils.h
        src/ParallelUtils.h
        src/PointDF.h
        src/SqlUtils.h
        src/StringUtils.h
        )

SET(${PROJECT_NAME}_LIB_SOURCES
        ${SOURCE_FILES}
        ${HEADER_FILES}
        src/PointDF.cpp
        src/PointDF.h
)

SET(DLL_LIB_INCLUDE_DIRS
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

INCLUDE_DIRECTORIES(
        DLL_LIB_INCLUDE_DIRS
)

FIND_PACKAGE(Qt5 COMPONENTS
        Core
        Test
        REQUIRED
        )

ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_LIB_SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
        Qt5::Core
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
