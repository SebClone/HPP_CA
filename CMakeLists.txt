cmake_minimum_required(VERSION 3.16)
project(encryptor LANGUAGES CXX)

# 1) C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 2) MPI finden
find_package(MPI REQUIRED)

# 3) Optional: OpenMP (für #pragma omp parallel/simd im Hauptloop)
find_package(OpenMP)

# 4) Globale Warnungen (optional, aber empfehlenswert)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif (MSVC)
    add_compile_options(/W4)
endif()

# 5) Include-Pfad für unsere Header
#    (erwarte: include/... enthält *.hpp inkl. hpp_rules_fast.hpp)
include_directories(${CMAKE_SOURCE_DIR}/include)

# 6) Kleine statische Libs für die Module
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

# 7) Neue schnelle Regeln (Template-Definition + explizite Instanziierungen)
add_library(hpp_rules_fast STATIC
    src/hpp_rules_fast.cpp
)
# hpp_rules_fast ruft collision()/reflection() aus hpp_rules auf → dagegen linken
target_link_libraries(hpp_rules_fast PUBLIC MPI::MPI_CXX hpp_rules)
target_include_directories(hpp_rules_fast PUBLIC ${CMAKE_SOURCE_DIR}/include)

# 8) Executable
add_executable(encryptor
    src/main.cpp
)
target_link_libraries(encryptor
    PRIVATE
        app_config
        io
        utilities
        hpp_rules
        hpp_rules_fast
        MPI::MPI_CXX
)

target_include_directories(encryptor PRIVATE ${CMAKE_SOURCE_DIR}/include)

# 9) OpenMP nur fürs Executable (dort sind die omp pragmas im Loop)
if (OpenMP_CXX_FOUND)
    target_link_libraries(encryptor PRIVATE OpenMP::OpenMP_CXX)
endif()

# 10) Release-Optimierungen / CPU-Vektorisierung
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    # -O3 + -march=native für Hot-Code (main & fast-rules)
    target_compile_options(encryptor     PRIVATE $<$<CONFIG:Release>:-O3 -march=native>)
    target_compile_options(hpp_rules_fast PRIVATE $<$<CONFIG:Release>:-O3 -march=native>)
elseif (MSVC)
    target_compile_options(encryptor      PRIVATE $<$<CONFIG:Release>:/O2>)
    target_compile_options(hpp_rules_fast PRIVATE $<$<CONFIG:Release>:/O2>)
    # Optional (falls deine Zielknoten AVX2 sicher unterstützen):
    # target_compile_options(encryptor      PRIVATE $<$<CONFIG:Release>:/arch:AVX2>)
    # target_compile_options(hpp_rules_fast PRIVATE $<$<CONFIG:Release>:/arch:AVX2>)
endif()

# 11) (optional) Sanitizer toggles:
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
if (ENABLE_ASAN AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(encryptor PRIVATE -fsanitize=address,undefined)
    target_link_options(encryptor PRIVATE -fsanitize=address,undefined)
endif();