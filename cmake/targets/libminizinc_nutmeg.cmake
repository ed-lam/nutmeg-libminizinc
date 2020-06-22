### MiniZinc Nutmeg Solver Target

if(NUTMEG_FOUND AND USE_NUTMEG)

  ### Compile target for the Nutmeg interface
  add_library(minizinc_nutmeg OBJECT
    solvers/nutmeg/nutmeg_constraints.cpp
    solvers/nutmeg/nutmeg_solverfactory.cpp
    solvers/nutmeg/nutmeg_solverinstance.cpp

    include/minizinc/solvers/nutmeg/nutmeg_constraints.hh
    include/minizinc/solvers/nutmeg_solverfactory.hh
    include/minizinc/solvers/nutmeg_solverinstance.hh
  )
  target_include_directories(minizinc_nutmeg PRIVATE "${NUTMEG_INCLUDE_DIRS}")
  #add_dependencies(minizinc_nutmeg minizinc_parser minizinc_geas minizinc_scip)
  add_dependencies(minizinc_nutmeg minizinc_parser)

  ### Setup correct compilation into the MiniZinc library
  target_compile_definitions(mzn PRIVATE HAS_NUTMEG)
  target_sources(mzn PRIVATE $<TARGET_OBJECTS:minizinc_nutmeg>)
  target_link_libraries(mzn ${NUTMEG_LIBRARIES})

endif()