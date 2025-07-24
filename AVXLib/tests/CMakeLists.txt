#https://github.com/KDAB/Charm/blob/master/Tests/CMakeLists.txt for example

SET(TEST_MODULE AVXLib)

SET(TEST_LIBRARIES
        ${TEST_MODULE}
        Qt5::Core
        Qt5::Test
    )


INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/bin)
MESSAGE(${CMAKE_BINARY_DIR}/bin)


if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm|arm64)")
    message(STATUS "Enabling NEON support via AvxToNeon")
    include_directories("${CMAKE_SOURCE_DIR}/ThirdPartyLibs/AvxToNeon")
    add_definitions(-DAVX2NEON)
    target_compile_options(${PROJECT_NAME} PRIVATE -march=armv8-a+fp+simd+crc)
else()
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
endif()

#Add classes to test here.
SET(TEST_CLASSES
        AVXUtilsTests
        )

foreach(class ${TEST_CLASSES})

    SET( CLASS_SRCS ${class}.cpp)
    ADD_EXECUTABLE( ${class} ${CLASS_SRCS} )
    TARGET_LINK_LIBRARIES( ${class} ${TEST_LIBRARIES} )
    SET_TARGET_PROPERTIES(${class} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
            ${CMAKE_BINARY_DIR}/bin )
    ADD_TEST( NAME ${class} COMMAND ${class} )

endforeach()
