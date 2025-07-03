cmake_minimum_required(VERSION 3.18)
PROJECT(EigenLib
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
        src/EigenKernelUtils.cpp
        src/EigenSparseUtils.cpp
        src/EigenUtils.cpp
        )

SET(HEADER_FILES
        src/EigenKernelUtils.h
        src/EigenSparseUtils.h
        src/EigenUtils.h
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
)

FIND_PACKAGE(Qt5 REQUIRED
        Test
	Network
)

#include(CheckCXXCompilerFlag)
## Check for AVX
#check_cxx_compiler_flag("-mavx" COMPILER_SUPPORTS_AVX)
#check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
#check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)
#
#if(COMPILER_SUPPORTS_AVX512)
#    message(STATUS "Enabling AVX-512 support, AVX, and AVX2 ${COMPILER_SUPPORTS_AVX512}")
#    add_compile_options(-mavx2 -mfma)
#    add_definitions(-DEIGEN_VECTORIZE_AVX)
#elseif(COMPILER_SUPPORTS_AVX2)
#    message(STATUS "Enabling AVX2 support")
#    add_definitions(-DEIGEN_VECTORIZE_AVX)
#    add_compile_options(-mavx2 -mfma)
#elseif(COMPILER_SUPPORTS_AVX)
#    message(STATUS "Enabling AVX support")
#    add_definitions(-DEIGEN_VECTORIZE_AVX)
#    add_compile_options(-mavx -mfma)
#elseif(MSVC)
#    message(STATUS "Enabling AVX support for MSVC")
#    add_compile_options(/arch:AVX) # Enable AVX support for MSVC
#    add_definitions(-DEIGEN_VECTORIZE_AVX)
#else()
#    message(WARNING "Your compiler does not support AVX, AVX2, or AVX-512 extensions.")
#endif()

ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_LIB_SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
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

string(TOUPPER ${PROJECT_NAME}_EXPORTS EXPORT_NAME)
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

ADD_SUBDIRECTORY("tests")
