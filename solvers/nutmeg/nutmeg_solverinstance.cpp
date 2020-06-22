/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Edward Lam <edward.lam@monash.edu>
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <minizinc/solvers/nutmeg_solverinstance.hh>
#include <minizinc/solvers/nutmeg/nutmeg_constraints.hh>

#include <minizinc/solvers/nutmeg_solverfactory.hh>
#include <minizinc/solvers/nutmeg_solverinstance.hh>

#include "scip/scip.h"

// ---------------------------------------------------------------------------------------------------------------------

using namespace Nutmeg;

namespace MiniZinc{
  NutmegSolverInstance::NutmegSolverInstance(Env &env, std::ostream &log, SolverInstanceBase::Options *opt)
      : SolverInstanceImpl<NutmegTypes>(env, log, opt), _flat(env.flat()),
        _solver(Nutmeg::Method::BC) {
    registerConstraints();
  }

  void NutmegSolverInstance::registerConstraint(std::string name, poster p) {
    _constraintRegistry.add("nutmeg_" + name, p);
    _constraintRegistry.add(name, p);
  }

  void NutmegSolverInstance::registerConstraints() {
    GCLock lock;

    /* Integer and Boolean Linear Constraints */
    registerConstraint("nutmeg_cp_int_bool_lin_le", NutmegConstraints::nutmeg_cp_int_bool_lin_le);
    registerConstraint("nutmeg_mip_int_bool_lin_le", NutmegConstraints::nutmeg_mip_int_bool_lin_le);
    registerConstraint("nutmeg_cp_int_bool_lin_eq", NutmegConstraints::nutmeg_cp_int_bool_lin_eq);
    registerConstraint("nutmeg_mip_int_bool_lin_eq", NutmegConstraints::nutmeg_mip_int_bool_lin_eq);

    /* Integer Linear Constraints */
    registerConstraint("nutmeg_cp_int_lin_le", NutmegConstraints::nutmeg_cp_int_lin_le);
    registerConstraint("nutmeg_mip_int_lin_le", NutmegConstraints::nutmeg_mip_int_lin_le);
    registerConstraint("nutmeg_cp_int_lin_le_imp", NutmegConstraints::nutmeg_cp_int_lin_le_imp);
    registerConstraint("nutmeg_mip_int_lin_le_imp", NutmegConstraints::nutmeg_mip_int_lin_le_imp);
    registerConstraint("nutmeg_cp_int_lin_le_reif", NutmegConstraints::nutmeg_cp_int_lin_le_reif);
    registerConstraint("nutmeg_mip_int_lin_le_reif", NutmegConstraints::nutmeg_mip_int_lin_le_reif);
    registerConstraint("nutmeg_cp_int_lin_eq", NutmegConstraints::nutmeg_cp_int_lin_eq);
    registerConstraint("nutmeg_mip_int_lin_eq", NutmegConstraints::nutmeg_mip_int_lin_eq);
    registerConstraint("nutmeg_cp_int_lin_eq_imp", NutmegConstraints::nutmeg_cp_int_lin_eq_imp);
    registerConstraint("nutmeg_mip_int_lin_eq_imp", NutmegConstraints::nutmeg_mip_int_lin_eq_imp);
    registerConstraint("nutmeg_cp_int_lin_eq_reif", NutmegConstraints::nutmeg_cp_int_lin_eq_reif);
    registerConstraint("nutmeg_cp_int_lin_ne", NutmegConstraints::nutmeg_cp_int_lin_ne);
    registerConstraint("nutmeg_cp_int_lin_ne_imp", NutmegConstraints::nutmeg_cp_int_lin_ne_imp);
    registerConstraint("nutmeg_cp_int_lin_ne_reif", NutmegConstraints::nutmeg_cp_int_lin_ne_reif);

    /* Integer Comparison Constraints */
    registerConstraint("nutmeg_cp_int_le", NutmegConstraints::nutmeg_cp_int_le);
    registerConstraint("nutmeg_cp_int_le_imp", NutmegConstraints::nutmeg_cp_int_le_imp);
    registerConstraint("nutmeg_cp_int_le_reif", NutmegConstraints::nutmeg_cp_int_le_reif);
    registerConstraint("int_eq", NutmegConstraints::int_eq); // SPECIAL
    registerConstraint("nutmeg_cp_int_eq", NutmegConstraints::nutmeg_cp_int_eq);
    registerConstraint("nutmeg_mip_int_eq", NutmegConstraints::nutmeg_mip_int_eq);
    registerConstraint("nutmeg_cp_int_eq_imp", NutmegConstraints::nutmeg_cp_int_eq_imp);
    registerConstraint("int_eq_reif", NutmegConstraints::int_eq_reif); // SPECIAL
    registerConstraint("nutmeg_cp_int_eq_reif", NutmegConstraints::nutmeg_cp_int_eq_reif);
    registerConstraint("nutmeg_cp_int_ne", NutmegConstraints::nutmeg_cp_int_ne);
    registerConstraint("nutmeg_cp_int_ne_imp", NutmegConstraints::nutmeg_cp_int_ne_imp);
    registerConstraint("nutmeg_cp_int_ne_reif", NutmegConstraints::nutmeg_cp_int_ne_reif);

    /* Integer Arithmetic Constraints */
    registerConstraint("nutmeg_cp_int_abs", NutmegConstraints::nutmeg_cp_int_abs);
    registerConstraint("nutmeg_cp_int_times", NutmegConstraints::nutmeg_cp_int_times);
    registerConstraint("nutmeg_cp_int_div", NutmegConstraints::nutmeg_cp_int_div);
    registerConstraint("nutmeg_cp_int_min", NutmegConstraints::nutmeg_cp_int_min);
    registerConstraint("nutmeg_cp_int_max", NutmegConstraints::nutmeg_cp_int_max);

    /* Boolean Linear Constraints */
    registerConstraint("nutmeg_cp_bool_lin_le", NutmegConstraints::nutmeg_cp_bool_lin_le);
    registerConstraint("nutmeg_mip_bool_lin_le", NutmegConstraints::nutmeg_mip_bool_lin_le);
    registerConstraint("nutmeg_cp_bool_lin_le_imp", NutmegConstraints::nutmeg_cp_bool_lin_le_imp);
    registerConstraint("nutmeg_mip_bool_lin_le_imp", NutmegConstraints::nutmeg_mip_bool_lin_le_imp);
    registerConstraint("nutmeg_cp_bool_lin_le_reif", NutmegConstraints::nutmeg_cp_bool_lin_le_reif);
    registerConstraint("nutmeg_mip_bool_lin_le_reif", NutmegConstraints::nutmeg_mip_bool_lin_le_reif);
    registerConstraint("nutmeg_cp_bool_lin_eq", NutmegConstraints::nutmeg_cp_bool_lin_eq);
    registerConstraint("nutmeg_mip_bool_lin_eq", NutmegConstraints::nutmeg_mip_bool_lin_eq);
    registerConstraint("nutmeg_cp_bool_lin_eq_imp", NutmegConstraints::nutmeg_cp_bool_lin_eq_imp);
    registerConstraint("nutmeg_mip_bool_lin_eq_imp", NutmegConstraints::nutmeg_mip_bool_lin_eq_imp);
    registerConstraint("nutmeg_cp_bool_lin_eq_reif", NutmegConstraints::nutmeg_cp_bool_lin_eq_reif);
    registerConstraint("nutmeg_cp_bool_lin_ne", NutmegConstraints::nutmeg_cp_bool_lin_ne);
    registerConstraint("nutmeg_cp_bool_lin_ne_imp", NutmegConstraints::nutmeg_cp_bool_lin_ne_imp);
    registerConstraint("nutmeg_cp_bool_lin_ne_reif", NutmegConstraints::nutmeg_cp_bool_lin_ne_reif);

    /* Boolean Arithmetic Constraints */
//    registerConstraint("nutmeg_cp_bool_or_imp", NutmegConstraints::nutmeg_cp_bool_or_imp);
//    registerConstraint("nutmeg_cp_bool_and_imp", NutmegConstraints::nutmeg_cp_bool_and_imp);
//    registerConstraint("nutmeg_cp_bool_xor_imp", NutmegConstraints::nutmeg_cp_bool_xor_imp);

    /* Clause Constraints */
    registerConstraint("nutmeg_cp_bool_clause", NutmegConstraints::nutmeg_cp_bool_clause);
//    registerConstraint("nutmeg_cp_bool_clause_imp", NutmegConstraints::nutmeg_cp_bool_clause_imp);
//    registerConstraint("nutmeg_cp_array_bool_or_imp", NutmegConstraints::nutmeg_cp_array_bool_or_imp);
//    registerConstraint("nutmeg_cp_array_bool_and_imp", NutmegConstraints::nutmeg_cp_array_bool_and_imp);
//    registerConstraint("nutmeg_cp_bool_clause_reif", NutmegConstraints::nutmeg_cp_bool_clause_reif);

    /* Coercion Constraints */
    registerConstraint("nutmeg_cp_bool2int", NutmegConstraints::nutmeg_cp_bool2int);
    registerConstraint("nutmeg_mip_bool2int", NutmegConstraints::nutmeg_mip_bool2int); // SPECIAL

    /* Element Constraints */
    registerConstraint("nutmeg_cp_array_int_element", NutmegConstraints::nutmeg_cp_array_int_element);
    registerConstraint("nutmeg_cp_array_var_int_element", NutmegConstraints::nutmeg_cp_array_var_int_element);
//    registerConstraint("nutmeg_cp_array_bool_element", NutmegConstraints::nutmeg_cp_array_bool_element);
//    registerConstraint("nutmeg_cp_array_var_bool_element", NutmegConstraints::nutmeg_cp_array_var_bool_element);

    /* Global Constraints */
    registerConstraint("nutmeg_cp_all_different", NutmegConstraints::nutmeg_cp_all_different);
    registerConstraint("nutmeg_cp_all_different_except_0", NutmegConstraints::nutmeg_cp_all_different_except_0);
//    registerConstraint("nutmeg_cp_at_most", NutmegConstraints::nutmeg_cp_at_most);
//    registerConstraint("nutmeg_cp_at_most1", NutmegConstraints::nutmeg_cp_at_most1);
    registerConstraint("nutmeg_cp_disjunctive", NutmegConstraints::nutmeg_cp_disjunctive);
    registerConstraint("nutmeg_cp_cumulative", NutmegConstraints::nutmeg_cp_cumulative);
    registerConstraint("nutmeg_cp_cumulative_var", NutmegConstraints::nutmeg_cp_cumulative_var);
    registerConstraint("nutmeg_cp_cumulative_optional", NutmegConstraints::nutmeg_cp_cumulative_optional);
    registerConstraint("nutmeg_cp_global_cardinality", NutmegConstraints::nutmeg_cp_global_cardinality);
    registerConstraint("nutmeg_cp_table_int", NutmegConstraints::nutmeg_cp_table_int);
  }

  void NutmegSolverInstance::processFlatZinc() {

      // Get options.
      auto _opt = static_cast<NutmegOptions&>(*_options);

      // Only support Boolean and integer variables.
      for (auto it = _flat->begin_vardecls(); it != _flat->end_vardecls(); ++it)
          if (!it->removed() && it->e()->type().isvar() && it->e()->type().dim() == 0)
          {
              VarDecl* vd = it->e();
              if (!vd->type().isint() && !vd->type().isbool())
              {
                  std::stringstream ssm;
                  ssm << "Type " << *vd->ti() << " is currently not supported by Nutmeg.";
                  throw InternalError(ssm.str());
              }
          }

      // Create integer variables.
      for (auto it = _flat->begin_vardecls(); it != _flat->end_vardecls(); ++it) {
          if (!it->removed() && it->e()->type().isvar() && it->e()->type().dim() == 0) {
              VarDecl* vd = it->e();
              if (vd->type().isint()) {
                  if (!vd->e()) {
                      Expression* domain = vd->ti()->domain();
                      if (domain) {
                          IntSetVal* isv = eval_intset(env().envi(), domain);

                          const auto min = static_cast<Int>(isv->min().toInt());
                          const auto max = static_cast<Int>(isv->max().toInt());
                          const auto& name = vd->id()->str().str();

                          IntVar var;
                          var = _solver.add_int_var(min, max, false, name);
                          if (isv->size() > 1) {
                              Vector<Int> vals;
                              for (int j = 0; j < isv->size(); ++j) {
                                  for (auto k = isv->min(j).toInt(); k <= isv->max(j).toInt(); ++k) {
                                      vals.push_back(static_cast<int>(k));
                                  }
                              }
                              debug_assert(!vals.empty());
                              _solver.add_indicator_vars(var, vals);
                          }
                          _variableMap.insert(vd->id(), NutmegVariable(var));
                      } else {
                          throw Error("NutmegSolverInstance::processFlatZinc: Error: Unbounded variable: " + vd->id()->str().str());
                      }
                  } else {
                      Expression* init = vd->e();
                      if (init->isa<Id>() || init->isa<ArrayAccess>()) {
                          NutmegVariable& var = resolveVar(init);
                          assert(var.isInt());
                          _variableMap.insert(vd->id(), NutmegVariable(var.intVar()));
                      } else {
                          auto il = init->cast<IntLit>()->v().toInt();

                          const auto constant = static_cast<Nutmeg::Int>(il);
                          auto var = _solver.add_int_var(constant, constant, false, fmt::format("constant_{}", constant));

                          _variableMap.insert(vd->id(), NutmegVariable(var));
                      }
                  }
              }
          }
      }

      // Create Boolean variables by substituting equality literals when encountering a reified equality constraint.
      for (ConstraintIterator it = _flat->begin_constraints(); it != _flat->end_constraints(); ++it)
          if (!it->removed())
              if (auto c = it->e()->dyn_cast<Call>())
              {
                  auto name = c->id().str();
                  if (name == "int_eq_reif")
                  {
                      auto x1 = c->arg(0);
                      auto x2 = c->arg(1);
                      auto x3 = c->arg(2);

                      if (x2->type().isvarint() && x1->type().ispar() && x1->type().isint())
                      {
                          std::swap(x1, x2);
                          debug_assert(x1->type().isvarint() && x2->type().ispar() && x2->type().isint());
                      }

                      if (x1->type().isvarint() && x2->type().ispar() && x2->type().isint())
                      {
                          // Find the Boolean variable declaration.
                          VarDecl* vd = nullptr;
                          for (auto it2 = _flat->begin_vardecls(); it2 != _flat->end_vardecls(); ++it2)
                              if (!it2->removed() && it2->e()->type().isvar() && it2->e()->type().dim() == 0)
                              {
                                  vd = it2->e();
                                  if (vd->id() == x3)
                                  {
                                      debug_assert(vd->type().isbool());
                                      break;
                                  }
                                  else
                                  {
                                      vd = nullptr;
                                  }
                              }

                          // Extract the literal.
                          if (vd)
                          {
                              // Get the variables.
                              auto int_var = asIntVar(x1);
                              const auto val = asInt(x2);
//                              std::cout << *x1 << " " << *x2 << " " << *x3 << std::endl;

                              // Proceed if the value is within the domain.
                              if (_solver.lb(int_var) <= val && val <= _solver.ub(int_var))
                              {
                                  // Create indicator variables.
                                  auto indicator_vars = _solver.add_indicator_vars(int_var);
                                  const auto lb = _solver.lb(int_var);
                                  const auto idx = val - lb;
                                  debug_assert(0 <= idx && idx < indicator_vars.size());
                                  auto lit = indicator_vars[idx];

                                  // Substitute the literal for the Boolean variable. If a substitution occurred for a
                                  // variable previously, add an equality constraint to equate the two literals.
                                  auto it2 = _variableMap.find(vd->id());
                                  if (it2 == _variableMap.end()) // Not substituted previously
                                  {
                                      _variableMap.insert(vd->id(), NutmegVariable(lit));

//                                      println("Substituting for {}", _solver.name(lit));
                                  }
                                  else if (!it2->second.boolVar().is_same(lit)) // Substituted previously and different
                                  {
                                      const auto new_lit = it2->second.boolVar();

                                      // Add CP equality constraint.
                                      {
                                          auto cp_lit = _solver.cp_var(lit);
                                          auto cp_new_lit = _solver.cp_var(new_lit);
                                          {
                                              vec<geas::clause_elt> clause;
                                              clause.push(cp_lit);
                                              clause.push(~cp_new_lit);
                                              if (!geas::add_clause(*_solver.cp_data(), clause))
                                              {
                                                  _solver.mark_as_infeasible();
                                              }
                                          }
                                          {
                                              vec<geas::clause_elt> clause;
                                              clause.push(~cp_lit);
                                              clause.push(cp_new_lit);
                                              if (!geas::add_clause(*_solver.cp_data(), clause))
                                              {
                                                  _solver.mark_as_infeasible();
                                              }
                                          }
                                      }

                                      // Add MIP equality constraint.
                                      {
                                          auto mip_lit = _solver.mip_var(lit);
                                          auto mip_new_lit = _solver.mip_var(new_lit);

                                          Vector<Float> mip_coeffs{1,-1};
                                          Vector<SCIP_VAR*> mip_vars{mip_lit, mip_new_lit};
                                          addMIPLinear(mip_coeffs, mip_vars, 0, 0);
                                      }

//                                      println("Equating literals {} and {}",
//                                              _solver.name(it2->second.boolVar()), _solver.name(lit));
                                  }

                                  // Reduce domain of the literal.
                                  if (!vd->e())
                                  {
                                      Expression* domain = vd->ti()->domain();
                                      if (domain)
                                      {
                                          IntBounds ib = compute_int_bounds(_env.envi(), domain);
                                          const auto lb = ib.l.toInt();
                                          const auto ub = ib.u.toInt();
                                          if (lb == ub)
                                          {
                                              fixVariable(lit, lb == 1);
                                          }
                                      }
                                  }
                                  else
                                  {
                                      Expression* init = vd->e();
                                      if (init->isa<Id>() || init->isa<ArrayAccess>())
                                      {
                                          not_yet_implemented();
                                      }
                                      else
                                      {
                                          const auto b = init->cast<BoolLit>()->v();
                                          fixVariable(lit, b);
                                      }
                                  }
                              }
//                              else
//                              {
//                                  std::cout << "Variable " << *x1 << " has domain " << _solver.lb(int_var) << "-" <<
//                                            _solver.ub(int_var) << " outside of value " << val << std::endl;
//                              }
                          }
                      }
                  }
              }

      // Fix some Boolean variables in reified equality constraints to false.
      for (ConstraintIterator it = _flat->begin_constraints(); it != _flat->end_constraints(); ++it)
          if (!it->removed())
              if (auto c = it->e()->dyn_cast<Call>())
              {
                  auto name = c->id().str();
                  if (name == "int_eq_reif")
                  {
                      auto x1 = c->arg(0);
                      auto x2 = c->arg(1);
                      auto x3 = c->arg(2);

                      if (x2->type().isvarint() && x1->type().ispar() && x1->type().isint())
                      {
                          std::swap(x1, x2);
                          debug_assert(x1->type().isvarint() && x2->type().ispar() && x2->type().isint());
                      }

                      if (x1->type().isvarint() && x2->type().ispar() && x2->type().isint())
                      {
                          // Find the Boolean variable declaration.
                          VarDecl* vd = nullptr;
                          for (auto it2 = _flat->begin_vardecls(); it2 != _flat->end_vardecls(); ++it2)
                              if (!it2->removed() && it2->e()->type().isvar() && it2->e()->type().dim() == 0)
                              {
                                  vd = it2->e();
                                  if (vd->id() == x3)
                                  {
                                      debug_assert(vd->type().isbool());
                                      break;
                                  }
                                  else
                                  {
                                      vd = nullptr;
                                  }
                              }

                          // Extract the literal.
                          if (vd)
                          {
                              // Get the variables.
                              auto int_var = asIntVar(x1);
                              const auto val = asInt(x2);
//                              std::cout << *x1 << " " << *x2 << " " << *x3 << std::endl;

                              // Proceed if the value is outside the domain.
                              if (!(_solver.lb(int_var) <= val && val <= _solver.ub(int_var)))
                              {
                                  auto it2 = _variableMap.find(vd->id());
                                  if (it2 == _variableMap.end())
                                  {
                                      // Substitute the false constant for the Boolean variable if a literal was not
                                      // added previously.
                                      _variableMap.insert(vd->id(), NutmegVariable(_solver.get_false()));
                                  }
                                  else
                                  {
                                      // Fix the literal to false if it was added previously.
                                      auto lit = it2->second.boolVar();
                                      fixVariable(lit, false);
                                  }
                              }
                          }
                      }
                  }
              }

      // Create the remaining Boolean variables.
      for (auto it = _flat->begin_vardecls(); it != _flat->end_vardecls(); ++it) {
          if (!it->removed() && it->e()->type().isvar() && it->e()->type().dim() == 0) {
              VarDecl* vd = it->e();
              if (vd->type().isbool()) {
                  if (_variableMap.find(vd->id()) == _variableMap.end()) {
                      if (!vd->e()) {
                          Expression* domain = vd->ti()->domain();
                          long long int lb, ub;
                          if (domain) {
                              IntBounds ib = compute_int_bounds(_env.envi(), domain);
                              lb = ib.l.toInt();
                              ub = ib.u.toInt();
                          } else {
                              lb = 0;
                              ub = 1;
                          }
                          if (lb == ub) {
                              BoolVar val = (lb == 0) ? _solver.get_false() : _solver.get_true();
                              _variableMap.insert(vd->id(), NutmegVariable(val));
                          } else {
                              const auto& name = vd->id()->str().str();
                              auto var = _solver.add_bool_var(name);
                              _variableMap.insert(vd->id(), NutmegVariable(var));
                          }
                      } else {
                          Expression* init = vd->e();
                          if (init->isa<Id>() || init->isa<ArrayAccess>()) {
                              NutmegVariable& var = resolveVar(init);
                              assert(var.isBool());
                              _variableMap.insert(vd->id(), NutmegVariable(var.boolVar()));
                          } else {
                              auto b = init->cast<BoolLit>()->v();
                              BoolVar val = b ? _solver.get_true() : _solver.get_false();
                              _variableMap.insert(vd->id(), NutmegVariable(val));
                          }
                      }
                  }
              }
          }
      }

      // Copy integer variables from Boolean variables in the MIP solver.
      for (ConstraintIterator it = _flat->begin_constraints(); it != _flat->end_constraints(); ++it) {
          if(!it->removed()) {
              if (auto c = it->e()->dyn_cast<Call>()) {
                  auto name = c->id().str();
                  if (name == "nutmeg_mip_bool2int") {
                      _constraintRegistry.post(c);
                  }
              }
          }
      }

      // Post constraints.
      for (ConstraintIterator it = _flat->begin_constraints(); it != _flat->end_constraints(); ++it) {
          if(!it->removed()) {
              if (auto c = it->e()->dyn_cast<Call>()) {
                  auto name = c->id().str();
                  if (name != "nutmeg_mip_bool2int") {
                      _constraintRegistry.post(c);
                  }
              }
          }
      }

      // Set objective.
      SolveI* si = _flat->solveItem();
      if(si->e()) {
          _obj_type = si->st();
          if (_obj_type == SolveI::ST_MIN) {
              _obj_var = std::unique_ptr<NutmegTypes::Variable>(new NutmegTypes::Variable(resolveVar(si->e())));
              _solver.add_mip_var(asIntVar(si->e()));
          } else if (_obj_type == SolveI::ST_MAX) {

              const auto original_var = asIntVar(si->e());
              _solver.add_mip_var(original_var);

              const auto var = _solver.add_int_var(-_solver.ub(original_var),
                                                   -_solver.lb(original_var),
                                                   true,
                                                   "-" + _solver.name(original_var));
              _obj_var = std::unique_ptr<NutmegTypes::Variable>(new NutmegTypes::Variable(var));

              release_assert(_solver.add_constr_linear({original_var, var}, {1, 1}, Sign::EQ, 0),
                             "Failed to create objective variable for maximization problem");
          }
      }
      if (!si->ann().isEmpty()) {
          fmt::print(stderr, "Warning: Nutmeg does not support search annotations\n");
      }
  }

  SolverInstanceBase::Status MiniZinc::NutmegSolverInstance::solve() {
      SolverInstanceBase::Status status = SolverInstance::ERROR;
      auto _opt = static_cast<NutmegOptions&>(*_options);

      // Solve if not failed at the root level.
      auto solver_status = _solver.get_status();
      if (solver_status != Nutmeg::Status::Infeasible)
      {
          // Get time limit.
          const auto time_limit = _opt.time == std::chrono::milliseconds(0) ?
                                  Nutmeg::Infinity :
                                  _opt.time.count() / 1e3;

          // Use MiniZinc output format.
          if (_opt.all_solutions)
          {
              if (_obj_type == SolveI::ST_SAT)
              {
                  fmt::print(stderr, "Warning: --all-solutions not yet implemented for satisfaction problems\n");
                  fflush(stderr);
              }
              else
              {
                  _solver.add_print_new_solution_function(std::bind(&NutmegSolverInstance::printSolution, this));
              }
          }

          // Solve.
          constexpr auto verbose = false;
          if (_obj_type == SolveI::ST_SAT)
          {
              _solver.satisfy(time_limit, verbose);
          }
          else
          {
              release_assert(_obj_var->isInt(), "Objective variable must be integer type");

              _solver.minimize(_obj_var->intVar(), time_limit, verbose);
          }
          solver_status = _solver.get_status();
      }

      // Process solution.
      switch (solver_status)
      {
          case Nutmeg::Status::Error:
              assert(false);
              status = SolverInstance::ERROR;
              break;
          case Nutmeg::Status::Infeasible:
              status = SolverInstance::UNSAT;
              break;
          case Nutmeg::Status::Unknown:
              status = SolverInstance::UNKNOWN;
              break;
          case Nutmeg::Status::Optimal:
              status = SolverInstance::OPT;
              if (_obj_type == SolveI::ST_SAT || !_opt.all_solutions)
              {
                  printSolution();
              }
              break;
          case Nutmeg::Status::Feasible:
              status = SolverInstance::SAT;
              if (_obj_type == SolveI::ST_SAT || !_opt.all_solutions)
              {
                  printSolution();
              }
              break;
      }

      // Print statistics for MiniZinc.
      if (_opt.statistics) {
          printStatistics(true);
      }

      // Exit.
      return status;
  }

  Expression* NutmegSolverInstance::getSolutionValue(Id* id) {
    id = id->decl()->id();
    if(id->type().isvar()) {
      NutmegVariable& var = resolveVar(id->decl()->id());
      switch (id->type().bt()) {
        case Type::BT_BOOL:
          assert(var.isBool());
          return constants().boollit(_solver.get_sol(var.boolVar()));
        case Type::BT_INT:
          assert(var.isInt());
          return IntLit::a(_solver.get_sol(var.intVar()));
        default:
          return nullptr;
      }
    } else {
      return id->decl()->e();
    }
  }

  void NutmegSolverInstance::resetSolver() {
      not_yet_implemented();
  }

  NutmegTypes::Variable& NutmegSolverInstance::resolveVar(Expression* e) {
    if (auto id = e->dyn_cast<Id>()) {
      return _variableMap.get(id->decl()->id());
    } else if (auto vd = e->dyn_cast<VarDecl>()) {
      return _variableMap.get(vd->id()->decl()->id());
    } else if (auto aa = e->dyn_cast<ArrayAccess>()) {
      auto ad = aa->v()->cast<Id>()->decl();
      auto idx = aa->idx()[0]->cast<IntLit>()->v().toInt();
      auto al = eval_array_lit(_env.envi(), ad->e());
      return _variableMap.get((*al)[idx]->cast<Id>());
    } else {
      std::stringstream ssm;
      ssm << "Expected Id, VarDecl or ArrayAccess instead of \"" << *e << "\"";
      throw InternalError(ssm.str());
    }
  }

  Nutmeg::Vector<bool> NutmegSolverInstance::asBool(ArrayLit* al) {
    Nutmeg::Vector<bool> vec(al->size());
    for (int i = 0; i < al->size(); ++i) {
      vec[i] = asBool((*al)[i]);
    }
    return vec;
  }

  Nutmeg::BoolVar NutmegSolverInstance::asBoolVar(Expression* e) {
    if (e->type().isvar()) {
      NutmegVariable& var = resolveVar(follow_id_to_decl(e));
      assert(var.isBool());
      return var.boolVar();
    } else {
      if(auto bl = e->dyn_cast<BoolLit>()) {
        return bl->v() ? _solver.get_true() : _solver.get_false();
      } else {
        std::stringstream ssm; ssm << "Expected bool or int literal instead of: " << *e;
        throw InternalError(ssm.str());
      }
    }
  }

  Nutmeg::Vector<Nutmeg::BoolVar> NutmegSolverInstance::asBoolVar(ArrayLit* al) {
    Nutmeg::Vector<Nutmeg::BoolVar> vec(al->size());
    for (int i = 0; i < al->size(); ++i) {
      vec[i] = this->asBoolVar((*al)[i]);
    }
    return vec;
  }

  Nutmeg::Vector<Nutmeg::Int> NutmegSolverInstance::asInt(ArrayLit* al) {
    Nutmeg::Vector<Nutmeg::Int> vec(al->size());
    for (int i = 0; i < al->size(); ++i) {
      vec[i] = this->asInt((*al)[i]);
    }
    return vec;
  }

  Nutmeg::IntVar NutmegSolverInstance::asIntVar(Expression* e) {
    if (e->type().isvar()) {
      NutmegVariable& var = resolveVar(follow_id_to_decl(e));
      assert(var.isInt());
      return var.intVar();
    } else {
      IntVal i;
      if(auto il = e->dyn_cast<IntLit>()) {
        i = il->v().toInt();
      } else if(auto bl = e->dyn_cast<BoolLit>()) {
        i = bl->v();
      } else {
        std::stringstream ssm; ssm << "Expected bool or int literal instead of: " << *e;
        throw InternalError(ssm.str());
      }
      if (i == 0) {
        return _solver.get_zero();
      } else {
        const auto constant = static_cast<Nutmeg::Int>(i.toInt());
        return _solver.add_int_var(constant, constant, false, fmt::format("constant_{}", constant));
      }
    }
  }

  Nutmeg::Vector<Nutmeg::IntVar> NutmegSolverInstance::asIntVar(ArrayLit* al) {
    Nutmeg::Vector<Nutmeg::IntVar> vec(al->size());
    for (int i = 0; i < al->size(); ++i) {
      vec[i] = this->asIntVar((*al)[i]);
    }
    return vec;
  }

  void NutmegSolverInstance::fixVariable(Nutmeg::BoolVar var, const bool val) {
      const auto lb = SCIPvarGetLbOriginal(_solver.mip_var(var));
      const auto ub = SCIPvarGetUbOriginal(_solver.mip_var(var));
//      println("Trying to fix {} to {} ({}-{})", _solver.name(var), val, lb, ub);
      if (lb > val || ub < val)
      {
          // Fixing variable to value outside bounds.
          _solver.mark_as_infeasible();
//          println("fixVariable making variable out of bounds (1)");
      }
      else
      {
          // Fix MIP vaiable.
          SCIP_Bool infeasible;
          SCIP_Bool fixed;
          scip_assert(SCIPfixVar(_solver.mip(), _solver.mip_var(var), val, &infeasible, &fixed));
//          println("fixVariable fixing variable {} to {}", _solver.name(var), val);
          if (infeasible)
          {
              _solver.mark_as_infeasible();
//              println("fixVariable making variable out of bounds (2)");
          }

          // Fix CP variable.
          if (!_solver.cp().post(val ? _solver.cp_var(var) : ~_solver.cp_var(var)))
          {
              _solver.mark_as_infeasible();
//              println("fixVariable making variable out of bounds (3)");
          }
      }
  }

  void NutmegSolverInstance::fixVariable(Nutmeg::IntVar var, const Nutmeg::Int val) {
      const auto lb = SCIPvarGetLbOriginal(_solver.mip_var(var));
      const auto ub = SCIPvarGetUbOriginal(_solver.mip_var(var));
//      println("Trying to fix {} to {} ({}-{})", _solver.name(var), val, lb, ub);
      if (lb > val || ub < val)
      {
          // Fixing variable to value outside bounds.
          _solver.mark_as_infeasible();
//          println("fixVariable making variable out of bounds (1)");
      }
      else
      {
          // Fix MIP vaiable.
          SCIP_Bool infeasible;
          SCIP_Bool fixed;
          scip_assert(SCIPfixVar(_solver.mip(), _solver.mip_var(var), val, &infeasible, &fixed));
//          println("fixVariable fixing variable {} to {}", _solver.name(var), val);
          if (infeasible)
          {
              _solver.mark_as_infeasible();
//              println("fixVariable making variable out of bounds (2)");
          }

          // Fix CP variable.
          if (!_solver.cp().post(_solver.cp_var(var) == val))
          {
              _solver.mark_as_infeasible();
//              println("fixVariable making variable out of bounds (3)");
          }
      }
  }

  void NutmegSolverInstance::addMIPLinear(Vector<Float> coeffs,
                                          Vector<SCIP_VAR*> vars,
                                          Float lhs,
                                          Float rhs)
  {
      // Check.
      release_assert(coeffs.size() == vars.size(), "Coefficients and variables arrays size do not match");
      release_assert(vars.size() > 0, "No variables to add to constraint");

      // Remove zero coefficients.
      for (Int idx = vars.size() - 1; idx >= 0; --idx)
          if (coeffs[idx] == 0.0)
          {
              vars.erase(vars.begin() + idx);
              coeffs.erase(coeffs.begin() + idx);
          }

      // Move fixed variables into the LHS and RHS.
      for (Int idx = vars.size() - 1; idx >= 0; --idx)
      {
          auto var = vars[idx];
          const auto lb = SCIPvarGetLbOriginal(var);
          const auto ub = SCIPvarGetUbOriginal(var);
          if (lb == ub)
          {
              const auto coeff = coeffs[idx];
              if (lhs > -INFINITY)
              {
                  lhs -= coeff * lb;
              }
              if (rhs < INFINITY)
              {
                  rhs -= coeff * lb;
              }
              vars.erase(vars.begin() + idx);
              coeffs.erase(coeffs.begin() + idx);
          }
      }

      // Change variable bounds if one variable. Otherwise add the constraint.
      if (vars.size() == 1 && coeffs[0] != 0.0)
      {
          auto var = vars[0];
          const auto coeff = coeffs[0];
          const auto lb = SCIPvarGetLbOriginal(var);
          const auto ub = SCIPvarGetUbOriginal(var);
          const auto new_lb = coeff > 0 ? lhs / coeff : rhs / coeff;
          const auto new_ub = coeff > 0 ? rhs / coeff : lhs / coeff;
          if (new_lb > lb)
          {
              scip_assert(SCIPchgVarLb(_solver.mip(), var, new_lb));
          }
          if (new_ub < ub)
          {
              scip_assert(SCIPchgVarUb(_solver.mip(), var, new_ub));
          }
      }
      else if (!(vars.empty() && lhs <= 0.0 && rhs >= 0.0))
      {
          SCIP_CONS* cons;
          scip_assert(SCIPcreateConsBasicLinear(_solver.mip(),
                                                &cons,
                                                fmt::format("linear_{}", _solver.nb_linear_constraints()++).c_str(),
                                                vars.size(),
                                                vars.data(),
                                                coeffs.data(),
                                                lhs,
                                                rhs));
          scip_assert(SCIPaddCons(_solver.mip(), cons));
          scip_assert(SCIPreleaseCons(_solver.mip(), &cons));
      }
  }

  void NutmegSolverInstance::addMIPIndicator(SCIP_VAR* r,
                                             Vector<Float>& coeffs,
                                             Vector<SCIP_VAR*>& vars,
                                             const Float lhs,
                                             const Float rhs)
  {
      // Check.
      release_assert(coeffs.size() == vars.size(), "Coefficients and variables arrays size do not match");
      release_assert(vars.size() > 0, "No variables to add to constraint");

      // Add linear constraint if r is already fixed.
      if (SCIPvarGetLbOriginal(r) == SCIPvarGetUbOriginal(r))
      {
          if (SCIPvarGetLbOriginal(r) == 1.0)
          {
              addMIPLinear(coeffs, vars, lhs, rhs);
          }
          else
          {
              release_assert(SCIPvarGetLbOriginal(r) == 0.0,
                             "Indicator variable {} is fixed to invalid value {}",
                             SCIPvarGetName(r), SCIPvarGetLbOriginal(r));
          }
      }
      else
      {
          // Add indicator constraint for <=.
          if (rhs < INFINITY)
          {
              SCIP_CONS* cons;
              scip_assert(SCIPcreateConsBasicIndicator(_solver.mip(),
                                                       &cons,
                                                       fmt::format("indicator_{}", _solver.nb_indicator_constraints()++).c_str(),
                                                       r,
                                                       vars.size(),
                                                       vars.data(),
                                                       coeffs.data(),
                                                       rhs));
              scip_assert(SCIPaddCons(_solver.mip(), cons));
              scip_assert(SCIPreleaseCons(_solver.mip(), &cons));
          }

          // Add indicator constraint for >=.
          if (lhs > -INFINITY)
          {
              auto neg_coeffs = coeffs;
              for (auto& x : neg_coeffs)
              {
                  x *= -1;
              }

              SCIP_CONS* cons;
              scip_assert(SCIPcreateConsBasicIndicator(_solver.mip(),
                                                       &cons,
                                                       fmt::format("indicator_{}", _solver.nb_indicator_constraints()++).c_str(),
                                                       r,
                                                       vars.size(),
                                                       vars.data(),
                                                       neg_coeffs.data(),
                                                       -lhs));
              scip_assert(SCIPaddCons(_solver.mip(), cons));
              scip_assert(SCIPreleaseCons(_solver.mip(), &cons));
          }
      }
  }

  void NutmegSolverInstance::printStatistics(bool fLegend) {
      auto& out = getSolns2Out()->getOutput();

      out << "%%%mzn-stat: solveTime=" << _solver.get_runtime() << std::endl;

      const auto solver_status = _solver.get_status();
      if (solver_status == Nutmeg::Status::Error)
      {
          out << "%%%mzn-stat: solveStatus=Error" << std::endl;
      }
      else if (solver_status == Nutmeg::Status::Infeasible)
      {
          out << "%%%mzn-stat: solveStatus=Infeasible" << std::endl;
      }
      else if (solver_status == Nutmeg::Status::Unknown)
      {
          out << "%%%mzn-stat: solveStatus=Unknown" << std::endl;
      }
      else if (solver_status == Nutmeg::Status::Optimal)
      {
          out << "%%%mzn-stat: solveStatus=Optimal" << std::endl;
      }
      else if (solver_status == Nutmeg::Status::Feasible)
      {
          out << "%%%mzn-stat: solveStatus=Feasible" << std::endl;
      }

      if (solver_status == Nutmeg::Status::Optimal || solver_status == Nutmeg::Status::Feasible)
      {
          auto dual_bound = _solver.get_dual_bound();
          auto primal_bound = _solver.get_primal_bound();
          if (_obj_type == SolveI::ST_MAX)
          {
              std::swap(primal_bound, dual_bound);
              dual_bound *= -1;
              primal_bound *= -1;
          }

          out << "%%%mzn-stat: objective=" << primal_bound << std::endl;
          out << "%%%mzn-stat: objectiveBound=" << dual_bound << std::endl;
      }

//      out << "%%%mzn-stat: nodes=" <<  << std::endl;
//      out << "%%%mzn-stat: solutions=" << st.solutions << std::endl;
  }

  Nutmeg_SolverFactory::Nutmeg_SolverFactory() {
    SolverConfig sc("org.minizinc.nutmeg", getVersion(nullptr));
    sc.name("Nutmeg");
    sc.mznlib("-Gnutmeg");
    sc.mznlibVersion(1);
    sc.supportsMzn(false);
    sc.description(getDescription(nullptr));
    sc.tags({"api","mip","cp","int","lcg",});
    sc.stdFlags({"-a", "-s", "-t"});
    sc.extraFlags({
    });
    SolverConfigs::registerBuiltinSolver(sc);
  };

  SolverInstanceBase::Options* Nutmeg_SolverFactory::createOptions() {
    return new NutmegOptions;
  }

  SolverInstanceBase* Nutmeg_SolverFactory::doCreateSI(Env& env, std::ostream& log, SolverInstanceBase::Options* opt) {
    return new NutmegSolverInstance(env, log, opt);
  }

  bool Nutmeg_SolverFactory::processOption(SolverInstanceBase::Options* opt, int &i, std::vector<std::string> &argv) {
    auto _opt = static_cast<NutmegOptions*>(opt);
    if (argv[i]=="-a" || argv[i]=="--all-solutions") {
      _opt->all_solutions = true;
//    } else if (argv[i]=="-f") {
//      _opt->free_search = true;
    } else if (argv[i]=="--solver-statistics" || argv[i]=="-s") {
      _opt->statistics = true;
    } else if (argv[i]=="--solver-time-limit" || argv[i]=="-t") {
      if (++i==argv.size()) return false;
      int time = atoi(argv[i].c_str());
      if(time >= 0)
        _opt->time = std::chrono::milliseconds(time);
    } else {
      return false;
    }
    return true;
  }

  void Nutmeg_SolverFactory::printHelp(std::ostream &os) {
    os << "Nutmeg solver plugin" << std::endl
       << std::endl;
  }
}

