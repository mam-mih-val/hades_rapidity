cmake_minimum_required(VERSION 3.11)
project(Rapidity)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
        ATTaskSkeleton
        GIT_REPOSITORY https://github.com/mam-mih-val/AnalysisTreeTaskSkeleton
        GIT_TAG        get_mass_check
)

FetchContent_MakeAvailable(ATTaskSkeleton)

find_package(ROOT REQUIRED COMPONENTS RIO Tree)

add_subdirectory(src)
