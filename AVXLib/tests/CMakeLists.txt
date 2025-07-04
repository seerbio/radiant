#https://github.com/KDAB/Charm/blob/master/Tests/CMakeLists.txt for example

SET(TEST_MODULE AVXLib)

SET(TEST_LIBRARIES
        ${TEST_MODULE}
        Qt5::Core
        Qt5::Test
    )


INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/bin)
MESSAGE(${CMAKE_BINARY_DIR}/bin)


if (COMPILER_SUPPORTS_AVX512)
    message(STATUS "Enabling AVX-512, AVX2, and AVX support")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2 -mfma")
elseif (COMPILER_SUPPORTS_AVX2)
    message(STATUS "Enabling AVX2 and AVX support")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2 -mfma")
elseif (COMPILER_SUPPORTS_AVX)
    message(STATUS "Enabling AVX support")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2 -mfma")
elseif (MSVC)
    message(STATUS "Enabling AVX support for MSVC")
#    target_compile_options(${PROJECT_NAME} PRIVATE /arch:AVX)
else()
    message(WARNING "AVX, AVX2, and AVX-512 are not supported by your compiler")
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
