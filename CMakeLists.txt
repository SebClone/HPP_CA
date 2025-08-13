cmake_minimum_required(VERSION 3.16)
project(encryptor LANGUAGES CXX)

# 1) C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 2) MPI finden
find_package(MPI REQUIRED)

# 3) Globale Warnungen (optional, aber empfehlenswert)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 4) Include-Pfad für unsere Header
include_directories(${CMAKE_SOURCE_DIR}/include)

# 5) Kleine statische Libs für die Module
add_library(app_config STATIC
    src/app_config.cpp
)
target_link_libraries(app_config PUBLIC MPI::MPI_CXX)
target_include_directories(app_config PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(io STATIC
    src/io.cpp
)
target_link_libraries(io PUBLIC MPI::MPI_CXX)
target_include_directories(io PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(utilities STATIC
    src/utilities.cpp
)
target_link_libraries(utilities PUBLIC MPI::MPI_CXX)
target_include_directories(utilities PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(hpp_rules STATIC
    src/hpp_rules.cpp
)
target_link_libraries(hpp_rules PUBLIC MPI::MPI_CXX)
target_include_directories(hpp_rules PUBLIC ${CMAKE_SOURCE_DIR}/include)

# 6) Executable (nur main.cpp)
add_executable(encryptor
    src/main.cpp
)
target_link_libraries(encryptor
    PRIVATE
        app_config
        io
        utilities
        hpp_rules
        MPI::MPI_CXX
)
target_include_directories(encryptor PRIVATE ${CMAKE_SOURCE_DIR}/include)

# (optional) Sanitizer toggles:
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
if (ENABLE_ASAN AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(encryptor PRIVATE -fsanitize=address,undefined)
    target_link_options(encryptor PRIVATE -fsanitize=address,undefined)
endif()
