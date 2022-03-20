cmake_minimum_required(VERSION 3.13)



###############################################################################
## PROJECT
## name and version
###############################################################################
project(Examples LANGUAGES CXX)



###############################################################################
## CONFIGURATION
## project configuration
###############################################################################
SET(${PARENT}_GENETICALGORITHM_ONEMAX_EXAMPLE           ON CACHE BOOL "" FORCE)
SET(${PARENT}_GENETICALGORITHM_CALIBTABLE_EXAMPLE       ON CACHE BOOL "" FORCE)



###############################################################################
## INCLUDING SUBDIRECTORIES
## Adding subdirectories according to the configuration
###############################################################################
if (${PARENT}_GENETICALGORITHM_ONEMAX_EXAMPLE)
    add_subdirectory(onemax)
endif()
if (${PARENT}_GENETICALGORITHM_ONEMAX_EXAMPLE)
    add_subdirectory(calibtable)
endif()
