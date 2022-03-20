cmake_minimum_required(VERSION 3.13)



###############################################################################
## PROJECT
## name and version
###############################################################################
project(GeneticAlgorithm VERSION 0.1.0 LANGUAGES CXX C)



###############################################################################
## SETTINGS
## basic project settings before use
###############################################################################
# If this project is used as a submodule, the variable should be overridden
# to "OFF" in the top-level application (to disable forced cache rewriting)
option(${PARENT}_SUBMODULE_CACHE_OVERWRITE "Enable forced cache rewriting" ON)
if (${PARENT}_SUBMODULE_CACHE_OVERWRITE)
    SET(REWRITE_FORCE "FORCE")
else()
    SET(REWRITE_FORCE "")
endif()



###############################################################################
## CONFIGURATION
## project configuration
###############################################################################
SET(${PARENT}_GENETICALGORITHM_EXAMMPLES     ON CACHE BOOL "" ${REWRITE_FORCE})
SET(${PARENT}_GENETICALGORITHM_APPLICATION   OFF CACHE BOOL "" ${REWRITE_FORCE})



###############################################################################
## INCLUDING SUBDIRECTORIES
## Adding subdirectories according to the configuration
###############################################################################
add_subdirectory(src)

if (${PARENT}_GENETICALGORITHM_EXAMMPLES)
    add_subdirectory(examples)
endif()
if (${PARENT}_GENETICALGORITHM_APPLICATION)
    add_subdirectory(applications)
endif()


