/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Edward Lam <edward.lam@monash.edu>
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <minizinc/solvers/nutmeg_solverfactory.hh>
#include <minizinc/solvers/nutmeg_solverinstance.hh>

namespace MiniZinc {
  namespace {
    void getWrapper() {
      static Nutmeg_SolverFactory _nutmeg_solverfactory;
    }
  }
  Nutmeg_SolverFactoryInitialiser::Nutmeg_SolverFactoryInitialiser() {
    getWrapper();
  }
}