/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Guido Tack <guido.tack@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <iostream>
#include <vector>

#include <minizinc/model.hh>
#include <minizinc/flatten.hh>
#include <minizinc/hash.hh>
#include <minizinc/statistics.hh>
#include <minizinc/solver_instance.hh>
#include <minizinc/solver_instance_defs.hh>
#include <minizinc/solns2out.hh>

namespace MiniZinc {

  /// An abstract SI
  class SolverInstanceBase {
  protected:
    Env& _env;
    Solns2Out* pS2Out=0;
    std::ostream& _log;
  public:
    
    /// Options base class
    /// A sub-class will be provided by each concrete SolverInstance
    class Options {
    public:
      bool verbose = false;
      bool printStatistics = false;
    };

    std::unique_ptr<Options> _options;
    
    typedef SolverInstance::Status Status;
    typedef SolverInstance::StatusReason StatusReason;
    Status _status;
    StatusReason _status_reason;
    
    SolverInstanceBase(Env& env, std::ostream& log, Options* options) : _env(env), _log(log), _options(options),
      _status(SolverInstance::UNKNOWN), _status_reason(SolverInstance::SR_OK) {}
    virtual ~SolverInstanceBase() { }
    
    /// Set/get the environment:
    virtual Env* getEnv() const { assert(&_env); return &_env; }
    virtual Env& env(void) const { return *getEnv(); }
    
    Solns2Out* getSolns2Out() const { assert(pS2Out); return pS2Out; }
    void setSolns2Out(Solns2Out* s2o) { pS2Out = s2o; }

    virtual void printSolution();
//     virtual void printSolution(ostream& );  // deprecated
    /// print statistics in form of comments
    virtual void printStatistics(bool fLegend=0) { }
    virtual void printStatisticsLine(bool fLegend=0) { }

    /// find the next solution
    virtual Status next(void) = 0;
    /// generate the solver-instance-representation from the flatzinc model
    virtual void processFlatZinc(void) = 0;
    /// clean up input model & flatzinc
    void cleanupForNonincrementalSolving() { getEnv()->envi().cleanupExceptOutput(); }
    /// solve the problem instance (according to the solve specification in the flatzinc model)
    virtual Status solve(void);
    /// return reason for status given by solve
    virtual StatusReason reason(void) {return _status_reason;}
    virtual Status status(void) {return _status;}
    

    /// reset the model to its core (removing temporary cts) and the solver to the root node of the search 
    void reset(void);
    /// reset the solver to the root node of the search TODO: difference between reset() and resetSolver()?
    virtual void resetSolver(void) = 0;
    /// reset the solver and add temporary constraints given by the iterator
    virtual void resetWithConstraints(Model::iterator begin, Model::iterator end);
    /// add permanent constraints given by the iterator to the solver instance
    virtual void processPermanentConstraints(Model::iterator begin, Model::iterator end);
  protected:
    /// flatten the search annotations, pushing them into the vector \a out
    void flattenSearchAnnotations(const Annotation& ann, std::vector<Expression*>& out);

    void flattenMultipleObjectives(const Annotation& ann, MultipleObjectives& mo);
    void flattenMultObjComponent(const Annotation& ann, MultipleObjectives::Objective& obj);

  private:
    SolverInstanceBase(const SolverInstanceBase&);
    SolverInstanceBase& operator= (const SolverInstanceBase&);
  };

  /// This implements a solver which is linked and returns its solution by assignSolutionToOutput()
  class SolverInstanceBase2 : public SolverInstanceBase {
  protected:
    virtual Expression* getSolutionValue(Id* id) = 0;

  public:
    /// Assign output for all vars: need public for callbacks
    // Default impl requires a Solns2Out object set up
    virtual void assignSolutionToOutput();
    /// Print solution to setup dest
    virtual void printSolution();
    
  protected:
    std::vector<VarDecl*> _varsWithOutput;    // this is to extract fzn vars. Identical to output()?  TODO

  public:
    SolverInstanceBase2(Env& env, std::ostream& log, SolverInstanceBase::Options* opt)
      : SolverInstanceBase(env, log, opt) {}
  };
  
  typedef void (*poster) (SolverInstanceBase&, const Call* call);
  class Registry {
  protected:
    std::unordered_map<std::string,poster> _registry;
    SolverInstanceBase& _base;
  public:
    Registry(SolverInstanceBase& base) : _base(base) {}
    void add(const std::string& name, poster p);
    void post(Call* c);      
    void cleanup() { _registry.clear(); }
  };

  /// Finally, this class also stores a mapping VarDecl->SolverVar and a constraint transformer
  /// It is a template holding parameterized VarId and Statistics, so cannot have members defined in a .cpp
  template<class Solver>
  class SolverInstanceImpl : public SolverInstanceBase2 {
  public:
    typename Solver::Statistics _statistics;
    virtual Statistics& getStatistics() { return _statistics; }

    typedef typename Solver::Variable VarId;

  protected:
    IdMap<VarId> _variableMap;      // this to find solver's variables given an Id
    Registry _constraintRegistry;

  public:
    SolverInstanceImpl(Env& env, std::ostream& log, SolverInstanceBase::Options* opt)
      : SolverInstanceBase2(env, log, opt), _constraintRegistry(*this) {}
    
  };

}
