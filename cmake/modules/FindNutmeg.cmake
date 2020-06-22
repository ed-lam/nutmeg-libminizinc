# - Try to find Nutmeg
# Once done this will define
#  NUTMEG_FOUND        - System has Nutmeg
#  NUTMEG_INCLUDE_DIRS - The Nutmeg include directories
#  NUTMEG_LIBRARIES    - The libraries needed to use Nutmeg
# User can set Nutmeg_ROOT to the preferred installation prefix
# Imported target Nutmeg will be created for linking purposes

# Find Nutmeg
# -----------

find_path(NUTMEG_INCLUDE Nutmeg/Nutmeg.h
          HINTS ${Nutmeg_ROOT} $ENV{Nutmeg_ROOT} ${NUTMEG_ROOT} $ENV{NUTMEG_ROOT} Nutmeg nutmeg
)
#message("NUTMEG NUTMEG_INCLUDE ${NUTMEG_INCLUDE}")

find_library(NUTMEG_LIBRARY NAMES nutmeg
             PATH_SUFFIXES lib build
             HINTS ${Nutmeg_ROOT} $ENV{Nutmeg_ROOT} ${NUTMEG_ROOT} $ENV{NUTMEG_ROOT} Nutmeg nutmeg
)
#message("NUTMEG NUTMEG_LIBRARY ${NUTMEG_LIBRARY}")

# Find SCIP
# ---------

if(NOT SCIP_INCLUDE)
  find_path(SCIP_INCLUDE scip/scip.h
            HINTS ${NUTMEG_INCLUDE}/scipoptsuite-6.0.2/scip/src
  )
#  message("NUTMEG SCIP_INCLUDE ${SCIP_INCLUDE}")
endif()

if(NOT SCIP_CONFIG_DIR)
  find_path(SCIP_CONFIG_DIR scip/config.h
            HINTS ${NUTMEG_INCLUDE}/build
  )
#  message("NUTMEG SCIP_CONFIG_DIR ${SCIP_CONFIG_DIR}")
endif()

if(NOT SCIP_LIBRARY)
  find_library(SCIP_LIBRARY NAMES scip
               HINTS ${NUTMEG_INCLUDE}/build/lib
  )
#  message("NUTMEG SCIP_LIBRARY ${SCIP_LIBRARY}")
endif()

# Find Geas
# ---------

if(NOT GEAS_INCLUDE)
  find_path(GEAS_INCLUDE geas/c/geas.h
            PATH_SUFFIXES include
            HINTS ${NUTMEG_INCLUDE}/geas
  )
#  message("NUTMEG GEAS_INCLUDE ${GEAS_INCLUDE}")
endif()

if(NOT GEAS_LIBRARY)
  find_library(GEAS_LIBRARY NAMES geas
               HINTS ${NUTMEG_INCLUDE}/build/geas
  )
#  message("NUTMEG GEAS_LIBRARY ${GEAS_LIBRARY}")
endif()

# Find CPLEX
# ----------

if(NOT CPLEX_INCLUDE)
  find_path(CPLEX_INCLUDE ilcplex/cplex.h
            HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR} ${CPLEX_DEFAULT_LOC}
            PATH_SUFFIXES include cplex/include)
  set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE})
#  message("NUTMEG CPLEX_INCLUDE ${CPLEX_INCLUDE}")
endif()

if(CPLEX_PLUGIN OR NOT CPLEX_LIB)
  set(CPLEX_LIBRARY "CPLEX_LIBRARY-NOTFOUND")
  foreach(CPLEX_LIB ${CPLEX_LIB_NAMES})
    find_library(CPLEX_LIBRARY NAMES cplex ${CPLEX_LIB}
                 HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR} ${CPLEX_DEFAULT_LOC}
                 PATH_SUFFIXES lib/x86-64_linux/static_pic lib/x86-64_osx/static_pic lib/x64_windows_vs2013/stat_mda cplex/lib/x86-64_linux/static_pic cplex/lib/x86-64_osx/static_pic cplex/lib/x64_windows_vs2013/stat_mda)
    if(NOT "${CPLEX_LIBRARY}" STREQUAL "CPLEX_LIBRARY-NOTFOUND")
      set(CPLEX_LIBRARIES ${CPLEX_LIBRARY})
#      message("NUTMEG CPLEX_LIBRARY ${CPLEX_LIBRARY}")
      break()
    endif()
  endforeach(CPLEX_LIB)
endif()

# Add Nutmeg
# ----------

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NUTMEG_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Nutmeg
  FOUND_VAR NUTMEG_FOUND
  REQUIRED_VARS NUTMEG_INCLUDE NUTMEG_LIBRARY
  FAIL_MESSAGE "Could NOT find Nutmeg, use Nutmeg_ROOT to hint its location"
)

mark_as_advanced(NUTMEG_INCLUDE NUTMEG_LIBRARY)

if(NUTMEG_FOUND)
  add_library(Nutmeg UNKNOWN IMPORTED)
  set_target_properties(Nutmeg PROPERTIES
    IMPORTED_LOCATION ${NUTMEG_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES "${NUTMEG_INCLUDE}"
  )
endif()

set(NUTMEG_LIBRARIES "${NUTMEG_LIBRARY}" "${SCIP_LIBRARY}" "${GEAS_LIBRARY}" "${CPLEX_LIBRARY}")
set(NUTMEG_INCLUDE_DIRS "${NUTMEG_INCLUDE}" "${NUTMEG_INCLUDE}/fmt/include" "${SCIP_INCLUDE}" "${SCIP_CONFIG_DIR}" "${GEAS_INCLUDE}" "${CPLEX_INCLUDE}")
