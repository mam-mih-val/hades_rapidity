cmake_minimum_required(VERSION 3.11)
project(Rapidity)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
        ATTaskSkeleton
        GIT_REPOSITORY https://github.com/eugene274/AnalysisTreeTaskSkeleton
        GIT_TAG        0.6.0
)

FetchContent_MakeAvailable(ATTaskSkeleton)

find_package(ROOT REQUIRED COMPONENTS RIO Tree)

add_subdirectory(src)
