/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Edward <edward.lam@monash.edu>
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <minizinc/flattener.hh>
#include <minizinc/solver.hh>
#include <Nutmeg/Nutmeg.h>

#include "scip/cons_linear.h"
#include "scip/cons_indicator.h"

namespace MiniZinc {

  class NutmegOptions : public SolverInstanceBase::Options {
  public:
    bool all_solutions = false;
    int conflicts = 0;
    bool free_search = false;
    int nr_solutions = 1;
    int obj_probe_limit = 0;
    bool statistics = false;
    std::chrono::milliseconds time = std::chrono::milliseconds(0);
  };

  class NutmegVariable {
  public:
    enum Type {BOOL_TYPE, INT_TYPE};
  protected:
    Type _t; // Type of the variable
    union {
      Nutmeg::BoolVar _bv;
      Nutmeg::IntVar _iv;
    };
  public:
    explicit NutmegVariable(const Nutmeg::BoolVar bv) : _t(BOOL_TYPE), _bv(bv) {};
    explicit NutmegVariable(const Nutmeg::IntVar iv) : _t(INT_TYPE), _iv(iv) {};

    NutmegVariable(const NutmegVariable& gv) : _t(gv._t) {
      switch (_t) {
        case BOOL_TYPE: _bv = gv._bv; break;
        case INT_TYPE: _iv = gv._iv; break;
      }
    }

    bool isBool() const { return _t == BOOL_TYPE; }
    bool isInt() const { return _t == INT_TYPE; }

    Nutmeg::BoolVar boolVar() { return _bv; }
    Nutmeg::IntVar intVar() { return _iv; }
  };

  class NutmegTypes {
  public:
    typedef NutmegVariable Variable;
    typedef MiniZinc::Statistics Statistics;
  };

  class NutmegSolverInstance : public SolverInstanceImpl<NutmegTypes> {
  public:
    NutmegSolverInstance(Env& env, std::ostream& log, SolverInstanceBase::Options* opt);
    ~NutmegSolverInstance() override = default;
    void processFlatZinc() override;
    Nutmeg::Model& solver() { return _solver; }

    void fixVariable(Nutmeg::BoolVar var, const bool val);
    void fixVariable(Nutmeg::IntVar var, const Nutmeg::Int val);
    void addMIPLinear(Nutmeg::Vector<Nutmeg::Float> coeffs,
                      Nutmeg::Vector<SCIP_VAR*> vars,
                      Nutmeg::Float lhs,
                      Nutmeg::Float rhs);
    void addMIPIndicator(SCIP_VAR* r,
                         Nutmeg::Vector<Nutmeg::Float>& coeffs,
                         Nutmeg::Vector<SCIP_VAR*>& vars,
                         const Nutmeg::Float lhs,
                         const Nutmeg::Float rhs);

    Status solve() override;
    Status next() override { return SolverInstance::ERROR; } // TODO: Implement
    void resetSolver() override;

    Expression* getSolutionValue(Id* id) override;
    void printStatistics(bool fLegend) override;

    // MiniZinc to Nutmeg conversions
    bool asBool(Expression* e) { return eval_bool(env().envi(), e); }
    Nutmeg::Vector<bool> asBool(ArrayLit* al);
    Nutmeg::BoolVar asBoolVar(Expression* e);
    Nutmeg::Vector<Nutmeg::BoolVar> asBoolVar(ArrayLit* al);
    Nutmeg::Vector<Nutmeg::Int> asInt(ArrayLit* al);
    Nutmeg::Int asInt(Expression* e) { return static_cast<Nutmeg::Int>(eval_int(env().envi(), e).toInt()); }
    Nutmeg::IntVar asIntVar(Expression* e);
    Nutmeg::Vector<Nutmeg::IntVar> asIntVar(ArrayLit* al);

  protected:
    Nutmeg::Model _solver;
    Model* _flat;

    SolveI::SolveType _obj_type = SolveI::ST_SAT;
    std::unique_ptr<NutmegTypes::Variable> _obj_var;

    NutmegTypes::Variable& resolveVar(Expression* e);

    void registerConstraint(std::string name, poster p);
    void registerConstraints();
  };

  class Nutmeg_SolverFactory: public SolverFactory {
  public:
    Nutmeg_SolverFactory();
    SolverInstanceBase::Options* createOptions() override;
    SolverInstanceBase* doCreateSI(Env& env, std::ostream& log, SolverInstanceBase::Options* opt) override;

    std::string getDescription(SolverInstanceBase::Options* opt) override { return "Nutmeg"; };
    std::string getVersion(SolverInstanceBase::Options* opt) override { return "0.0.2"; }
    std::string getId() override { return "org.minizinc.nutmeg"; }

    bool processOption(SolverInstanceBase::Options* opt, int& i, std::vector<std::string>& argv) override;
    void printHelp(std::ostream& os) override;
  };

}
