/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Edward Lam <edward.lam@monash.edu>
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-pro-type-static-cast-downcast"

#include <minizinc/solvers/nutmeg/nutmeg_constraints.hh>
#include <minizinc/solvers/nutmeg_solverinstance.hh>

#include <geas/constraints/builtins.h>
#include <geas/constraints/flow/flow.h>

using namespace Nutmeg;

namespace MiniZinc
{
namespace NutmegConstraints
{

#define SI static_cast<NutmegSolverInstance&>(s)
#define SOL SI.solver()
#define EXPR(X) call->arg(X)
#define BOOL(X) SI.asBool(EXPR(X))
#define BOOLARRAY(X) SI.asBool(ARRAY(X))
#define BOOLVAR(X) SI.asBoolVar(EXPR(X))
#define BOOLVARARRAY(X) SI.asBoolVar(ARRAY(X))
#define INT(X) SI.asInt(EXPR(X))
#define INTARRAY(X) SI.asInt(ARRAY(X))
#define INTVAR(X) SI.asIntVar(EXPR(X))
#define INTVARARRAY(X) SI.asIntVar(ARRAY(X))
#define PAR(X) call->arg(X)->type().ispar()
#define ARRAY(X) eval_array_lit(s.env().envi(), call->arg(X))

#define CP SOL.cp()
#define CPDATA SOL.cp_data()
#define MIP SOL.mip()
#define CPVAR(X) SOL.cp_var(X)
#define MIPVAR(X) SOL.mip_var(SOL.add_mip_var(X))
#define CPNEGVAR(X) ~CPVAR(X)
#define MIPNEGVAR(X) MIPVAR(SOL.get_neg(X))
#define ZEROVAR SOL.get_zero()
#define INFINITY SCIPinfinity(MIP)
#define CPCONSTRAINT(CONSTR, ...) { if (!geas::CONSTR(CPDATA, ##__VA_ARGS__)) SOL.mark_as_infeasible(); }
#define CPCLAUSE(...) { if (!geas::add_clause(*CPDATA, ##__VA_ARGS__)) SOL.mark_as_infeasible(); }

// Integer and Boolean Linear Constraints
// --------------------------------------

void nutmeg_cp_int_bool_lin_le(SolverInstanceBase& s, const Call* call)
{
    const auto& x = INTVAR(0);
    const auto& bool_coeffs = INTARRAY(1);
    const auto& bool_vars = BOOLVARARRAY(2);
    const auto rhs = INT(3);

    vec<int> cp_bool_coeffs;
    vec<int> cp_bool_neg_coeffs;
    vec<geas::patom_t> cp_bool_vars;
    for (size_t i = 0; i < bool_vars.size(); ++i)
    {
        cp_bool_coeffs.push(bool_coeffs[i]);
        cp_bool_neg_coeffs.push(-bool_coeffs[i]);
        cp_bool_vars.push(CPVAR(bool_vars[i]));
    }
    // true -> sum(coeffs[i] * vars[i]) + x <= rhs
    // true -> sum(coeffs[i] * vars[i]) - rhs <= -x
    // true -> -x >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, geas::at_True, -CPVAR(x), cp_bool_coeffs, cp_bool_vars, -rhs);
}
void nutmeg_mip_int_bool_lin_le(SolverInstanceBase& s, const Call* call)
{
    const auto& int_coeffs = INTARRAY(0);
    const auto& int_vars = INTVARARRAY(1);
    const auto& bool_coeffs = INTARRAY(2);
    const auto& bool_vars = BOOLVARARRAY(3);
    const auto rhs = INT(4);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < int_vars.size(); ++i)
    {
        mip_coeffs.push_back(int_coeffs[i]);
        mip_vars.push_back(MIPVAR(int_vars[i]));
    }
    for (size_t i = 0; i < bool_vars.size(); ++i)
    {
        mip_coeffs.push_back(bool_coeffs[i]);
        mip_vars.push_back(MIPVAR(bool_vars[i]));
    }
    SI.addMIPLinear(mip_coeffs, mip_vars, -INFINITY, rhs);
}

void nutmeg_cp_int_bool_lin_eq(SolverInstanceBase& s, const Call* call)
{
    const auto& x = INTVAR(0);
    const auto& bool_coeffs = INTARRAY(1);
    const auto& bool_vars = BOOLVARARRAY(2);
    const auto rhs = INT(3);

    vec<int> cp_bool_coeffs;
    vec<int> cp_bool_neg_coeffs;
    vec<geas::patom_t> cp_bool_vars;
    for (size_t i = 0; i < bool_vars.size(); ++i)
    {
        cp_bool_coeffs.push(bool_coeffs[i]);
        cp_bool_neg_coeffs.push(-bool_coeffs[i]);
        cp_bool_vars.push(CPVAR(bool_vars[i]));
    }
    // true -> sum(coeffs[i] * vars[i]) + x <= rhs
    // true -> sum(coeffs[i] * vars[i]) - rhs <= -x
    // true -> -x >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, geas::at_True, -CPVAR(x), cp_bool_coeffs, cp_bool_vars, -rhs);
    // true -> sum(coeffs[i] * vars[i]) + x >= rhs
    // true -> sum(coeffs[i] * vars[i]) - rhs >= -x
    // true -> sum(-coeffs[i] * vars[i]) + rhs <= x
    // true -> x >= sum(-coeffs[i] * vars[i]) + rhs
    CPCONSTRAINT(bool_linear_ge, geas::at_True, CPVAR(x), cp_bool_neg_coeffs, cp_bool_vars, rhs);
}
void nutmeg_mip_int_bool_lin_eq(SolverInstanceBase& s, const Call* call)
{
    const auto& int_coeffs = INTARRAY(0);
    const auto& int_vars = INTVARARRAY(1);
    const auto& bool_coeffs = INTARRAY(2);
    const auto& bool_vars = BOOLVARARRAY(3);
    const auto rhs = INT(4);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < int_vars.size(); ++i)
    {
        mip_coeffs.push_back(int_coeffs[i]);
        mip_vars.push_back(MIPVAR(int_vars[i]));
    }
    for (size_t i = 0; i < bool_vars.size(); ++i)
    {
        mip_coeffs.push_back(bool_coeffs[i]);
        mip_vars.push_back(MIPVAR(bool_vars[i]));
    }
    SI.addMIPLinear(mip_coeffs, mip_vars, rhs, rhs);
}

// Integer Linear Constraints
// --------------------------

void nutmeg_cp_int_lin_le(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);

    vec<int> cp_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_le, cp_coeffs, cp_vars, rhs);
}
void nutmeg_mip_int_lin_le(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPLinear(mip_coeffs, mip_vars, -INFINITY, rhs);
}

void nutmeg_cp_int_lin_le_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_le, cp_coeffs, cp_vars, rhs, CPVAR(r));
}
void nutmeg_mip_int_lin_le_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPIndicator(MIPVAR(r), mip_coeffs, mip_vars, -INFINITY, rhs);
}

void nutmeg_cp_int_lin_le_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_le, cp_coeffs, cp_vars, rhs, CPVAR(r));
    CPCONSTRAINT(linear_le, cp_neg_coeffs, cp_vars, -(rhs + 1), ~CPVAR(r));
}
void nutmeg_mip_int_lin_le_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPIndicator(MIPVAR(r), mip_coeffs, mip_vars, -INFINITY, rhs);
    SI.addMIPIndicator(MIPNEGVAR(r), mip_coeffs, mip_vars, rhs + 1, INFINITY);
}

void nutmeg_cp_int_lin_eq(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_le, cp_coeffs, cp_vars, rhs);
    CPCONSTRAINT(linear_le, cp_neg_coeffs, cp_vars, -rhs);
}
void nutmeg_mip_int_lin_eq(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPLinear(mip_coeffs, mip_vars, rhs, rhs);
}

void nutmeg_cp_int_lin_eq_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_le, cp_coeffs, cp_vars, rhs, CPVAR(r));
    CPCONSTRAINT(linear_le, cp_neg_coeffs, cp_vars, -rhs, CPVAR(r));
}
void nutmeg_mip_int_lin_eq_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPIndicator(MIPVAR(r), mip_coeffs, mip_vars, rhs, rhs);
}

void nutmeg_cp_int_lin_eq_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_le, cp_coeffs, cp_vars, rhs, CPVAR(r));
    CPCONSTRAINT(linear_le, cp_neg_coeffs, cp_vars, -rhs, CPVAR(r));
    CPCONSTRAINT(linear_ne, cp_coeffs, cp_vars, rhs, ~CPVAR(r));
}

void nutmeg_cp_int_lin_ne(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);

    vec<int> cp_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_ne, cp_coeffs, cp_vars, rhs);
}

void nutmeg_cp_int_lin_ne_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_ne, cp_coeffs, cp_vars, rhs, CPVAR(r));
}

void nutmeg_cp_int_lin_ne_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = INTVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(linear_ne, cp_coeffs, cp_vars, rhs, CPVAR(r));
    CPCONSTRAINT(linear_le, cp_coeffs, cp_vars, rhs, ~CPVAR(r));
    CPCONSTRAINT(linear_le, cp_neg_coeffs, cp_vars, -rhs, ~CPVAR(r));
}

// Integer Comparison Constraints
// ------------------------------

void nutmeg_cp_int_le(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_le, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), INT(2));
}

void nutmeg_cp_int_le_imp(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_le, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), INT(2), CPVAR(BOOLVAR(3)));
}

void nutmeg_cp_int_le_reif(SolverInstanceBase& s, const Call* call) {
    CPCONSTRAINT(int_le, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), INT(2), CPVAR(BOOLVAR(3)));
    CPCONSTRAINT(int_le, CPVAR(INTVAR(1)), CPVAR(INTVAR(0)), -(INT(2)+1), ~CPVAR(BOOLVAR(3)));
}

void int_eq(SolverInstanceBase& s, const Call* call) // SPECIAL
{
    nutmeg_cp_int_eq(s, call);
    nutmeg_mip_int_eq(s, call);
}
void nutmeg_cp_int_eq(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_eq, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)));
}
void nutmeg_mip_int_eq(SolverInstanceBase& s, const Call* call)
{
    Vector<Float> mip_coeffs{1,-1};
    Vector<SCIP_VAR*> mip_vars{MIPVAR(INTVAR(0)), MIPVAR(INTVAR(1))};
    SI.addMIPLinear(mip_coeffs, mip_vars, 0, 0);
}

void nutmeg_cp_int_eq_imp(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_eq, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), CPVAR(BOOLVAR(2)));
}

void int_eq_reif(SolverInstanceBase& s, const Call* call)
{
    // a == b <-> r
    if (PAR(0) && PAR(1)) // a and b are fixed
    {
        const auto a = INT(0);
        const auto b = INT(1);
        if (PAR(2)) // r is also fixed
        {
            const auto r = BOOL(2);
            if ((r && a != b) || (!r && a == b))
            {
                SOL.mark_as_infeasible();
            }
        }
        else // r is not fixed
        {
            const auto r = BOOLVAR(2);
//            println("int_eq_reif({},{},{})", a, b, SOL.name(r));
            SI.fixVariable(r, a == b);
        }
    }
    else if ((!PAR(0) && PAR(1)) || (PAR(0) && !PAR(1)))  // One of a or b is fixed
    {
        const auto var = !PAR(0) && PAR(1) ? INTVAR(0) : INTVAR(1);
        const auto val = !PAR(0) && PAR(1) ? INT(1) : INT(0);
        if (PAR(2)) // r is also fixed
        {
            const auto r = BOOL(2);
            if (r) // r is true
            {
                // Fix variable a to constant b.
                SI.fixVariable(var, val);
            }
            else
            {
                // Create CP reified constraint.
                nutmeg_cp_int_eq_reif(s, call);
            }
        }
        else // r is not fixed
        {
            const auto r = BOOLVAR(2);
//            println("int_eq_reif({},{},{}) {} {}", SOL.name(var), val, SOL.name(r), SOL.lb(var), SOL.ub(var));

            if (SOL.lb(var) == SOL.ub(var)) // a is fixed
            {
                SI.fixVariable(r, SOL.lb(var) == val);
            }
            else // a and r are not fixed
            {
                release_assert(
                    (SOL.lb(var) <= val && val <= SOL.ub(var) && SOL.has_mip_indicator_vars(var)) ||
                    ((val > SOL.ub(var) || val < SOL.lb(var)) && r.is_same(SOL.get_false())),
                    "Internal error in int_eq_reif constraint"
                );
            }
        }
    }
    else  // Neither x nor y are fixed
    {
//        println("int_eq_reif({},{},{})", SOL.name(INTVAR(0)), SOL.name(INTVAR(1)), SOL.name(BOOLVAR(2)));

        // Create CP reified constraint.
        nutmeg_cp_int_eq_reif(s, call);

        // Create MIP relaxation using half reification r -> x1 - x2 == 0
        Vector<Float> mip_coeffs{1, -1};
        Vector<SCIP_VAR*> mip_vars{MIPVAR(INTVAR(0)), MIPVAR(INTVAR(1))};
        SI.addMIPIndicator(MIPVAR(BOOLVAR(2)), mip_coeffs, mip_vars, 0, 0);
    }
}
void nutmeg_cp_int_eq_reif(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_eq, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), CPVAR(BOOLVAR(2)));
    CPCONSTRAINT(int_ne, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), ~CPVAR(BOOLVAR(2)));
}

void nutmeg_cp_int_ne(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_ne, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)));
}

void nutmeg_cp_int_ne_imp(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_ne, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), CPVAR(BOOLVAR(2)));
}

void nutmeg_cp_int_ne_reif(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_ne, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), CPVAR(BOOLVAR(2)));
    CPCONSTRAINT(int_eq, CPVAR(INTVAR(0)), CPVAR(INTVAR(1)), ~CPVAR(BOOLVAR(2)));
}

// Integer Arithmetic Constraints
// ------------------------------

void nutmeg_cp_int_abs(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_abs, CPVAR(INTVAR(1)), CPVAR(INTVAR(0)));
}

void nutmeg_cp_int_times(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_mul, CPVAR(INTVAR(2)), CPVAR(INTVAR(0)), CPVAR(INTVAR(1)));
}

void nutmeg_cp_int_div(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(int_div, CPVAR(INTVAR(2)), CPVAR(INTVAR(0)), CPVAR(INTVAR(1)));
}

void nutmeg_cp_int_min(SolverInstanceBase& s, const Call* call)
{
    vec<geas::intvar> vars = {-CPVAR(INTVAR(0)), -CPVAR(INTVAR(1))};
    CPCONSTRAINT(int_max, -CPVAR(INTVAR(2)), vars);
}

void nutmeg_cp_int_max(SolverInstanceBase& s, const Call* call)
{
    vec<geas::intvar> vars = {CPVAR(INTVAR(0)), CPVAR(INTVAR(1))};
    CPCONSTRAINT(int_max, CPVAR(INTVAR(2)), vars);
}

// Boolean Linear Constraints
// --------------------------

void nutmeg_cp_bool_lin_le(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);

    vec<int> cp_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    // true -> sum(coeffs[i] * vars[i]) <= rhs
    // true -> sum(coeffs[i] * vars[i]) - rhs <= 0
    // true -> 0 >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, geas::at_True, CPVAR(ZEROVAR), cp_coeffs, cp_vars, -rhs);
}
void nutmeg_mip_bool_lin_le(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPLinear(mip_coeffs, mip_vars, -INFINITY, rhs);
}

void nutmeg_cp_bool_lin_le_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    // r -> sum(coeffs[i] * vars[i]) <= rhs
    // r -> sum(coeffs[i] * vars[i]) - rhs <= 0
    // r -> 0 >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, CPVAR(r), CPVAR(ZEROVAR), cp_coeffs, cp_vars, -rhs);
}
void nutmeg_mip_bool_lin_le_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPIndicator(MIPVAR(r), mip_coeffs, mip_vars, -INFINITY, rhs);
}

void nutmeg_cp_bool_lin_le_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    // r -> sum(coeffs[i] * vars[i]) <= rhs
    // r -> sum(coeffs[i] * vars[i]) - rhs <= 0
    // r -> 0 >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, CPVAR(r), CPVAR(ZEROVAR), cp_coeffs, cp_vars, -rhs);
    // ~r -> sum(coeffs[i] * vars[i]) >= rhs + 1
    // ~r -> sum(coeffs[i] * vars[i]) - (rhs + 1) >= 0
    // ~r -> sum(-coeffs[i] * vars[i]) + (rhs + 1) <= 0
    // ~r -> 0 >= sum(-coeffs[i] * vars[i]) + (rhs + 1)
    CPCONSTRAINT(bool_linear_ge, ~CPVAR(r), CPVAR(ZEROVAR), cp_neg_coeffs, cp_vars, rhs + 1);
}
void nutmeg_mip_bool_lin_le_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPIndicator(MIPVAR(r), mip_coeffs, mip_vars, -INFINITY, rhs);
    SI.addMIPIndicator(MIPNEGVAR(r), mip_coeffs, mip_vars, rhs + 1, INFINITY);
}

void nutmeg_cp_bool_lin_eq(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    // true -> sum(coeffs[i] * vars[i]) <= rhs
    // true -> sum(coeffs[i] * vars[i]) - rhs <= 0
    // true -> 0 >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, geas::at_True, CPVAR(ZEROVAR), cp_coeffs, cp_vars, -rhs);
    // true -> sum(coeffs[i] * vars[i]) >= rhs
    // true -> sum(coeffs[i] * vars[i]) - rhs >= 0
    // true -> sum(-coeffs[i] * vars[i]) + rhs <= 0
    // true -> 0 >= sum(-coeffs[i] * vars[i]) + rhs
    CPCONSTRAINT(bool_linear_ge, geas::at_True, CPVAR(ZEROVAR), cp_neg_coeffs, cp_vars, rhs);
}
void nutmeg_mip_bool_lin_eq(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPLinear(mip_coeffs, mip_vars, rhs, rhs);
}

void nutmeg_cp_bool_lin_eq_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    // r -> sum(coeffs[i] * vars[i]) <= rhs
    // r -> sum(coeffs[i] * vars[i]) - rhs <= 0
    // r -> 0 >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, CPVAR(r), CPVAR(ZEROVAR), cp_coeffs, cp_vars, -rhs);
    // r -> sum(coeffs[i] * vars[i]) >= rhs
    // r -> sum(coeffs[i] * vars[i]) - rhs >= 0
    // r -> sum(-coeffs[i] * vars[i]) + rhs <= 0
    // r -> 0 >= sum(-coeffs[i] * vars[i]) + rhs
    CPCONSTRAINT(bool_linear_ge, CPVAR(r), CPVAR(ZEROVAR), cp_neg_coeffs, cp_vars, rhs);
}
void nutmeg_mip_bool_lin_eq_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    Vector<Float> mip_coeffs;
    Vector<SCIP_VAR*> mip_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        mip_coeffs.push_back(coeffs[i]);
        mip_vars.push_back(MIPVAR(vars[i]));
    }
    SI.addMIPIndicator(MIPVAR(r), mip_coeffs, mip_vars, rhs, rhs);
}

void nutmeg_cp_bool_lin_eq_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    // r -> sum(coeffs[i] * vars[i]) <= rhs
    // r -> sum(coeffs[i] * vars[i]) - rhs <= 0
    // r -> 0 >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, CPVAR(r), CPVAR(ZEROVAR), cp_coeffs, cp_vars, -rhs);
    // r -> sum(coeffs[i] * vars[i]) >= rhs
    // r -> sum(coeffs[i] * vars[i]) - rhs >= 0
    // r -> sum(-coeffs[i] * vars[i]) + rhs <= 0
    // r -> 0 >= sum(-coeffs[i] * vars[i]) + rhs
    CPCONSTRAINT(bool_linear_ge, CPVAR(r), CPVAR(ZEROVAR), cp_neg_coeffs, cp_vars, rhs);
    // ~r -> sum(coeffs[i] * vars[i]) != rhs
    not_yet_implemented();
    CPCONSTRAINT(bool_linear_ne, cp_coeffs, cp_vars, rhs, ~CPVAR(r));
}

void nutmeg_cp_bool_lin_ne(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);

    vec<int> cp_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    not_yet_implemented();
    CPCONSTRAINT(bool_linear_ne, cp_coeffs, cp_vars, rhs);
}

void nutmeg_cp_bool_lin_ne_imp(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    not_yet_implemented();
    CPCONSTRAINT(bool_linear_ne, cp_coeffs, cp_vars, rhs, CPVAR(r));
}

void nutmeg_cp_bool_lin_ne_reif(SolverInstanceBase& s, const Call* call)
{
    const auto& coeffs = INTARRAY(0);
    const auto& vars = BOOLVARARRAY(1);
    const auto rhs = INT(2);
    const auto& r = BOOLVAR(3);

    vec<int> cp_coeffs;
    vec<int> cp_neg_coeffs;
    vec<geas::patom_t> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_coeffs.push(coeffs[i]);
        cp_neg_coeffs.push(-coeffs[i]);
        cp_vars.push(CPVAR(vars[i]));
    }
    // r -> sum(coeffs[i] * vars[i]) != rhs
    not_yet_implemented();
    CPCONSTRAINT(bool_linear_ne, cp_coeffs, cp_vars, rhs, CPVAR(r));
    // ~r -> sum(coeffs[i] * vars[i]) <= rhs
    // ~r -> sum(coeffs[i] * vars[i]) - rhs <= 0
    // ~r -> 0 >= sum(coeffs[i] * vars[i]) - rhs
    CPCONSTRAINT(bool_linear_ge, ~CPVAR(r), CPVAR(ZEROVAR), cp_coeffs, cp_vars, -rhs);
    // ~r -> sum(coeffs[i] * vars[i]) >= rhs
    // ~r -> sum(coeffs[i] * vars[i]) - rhs >= 0
    // ~r -> sum(-coeffs[i] * vars[i]) + rhs <= 0
    // ~r -> 0 >= sum(-coeffs[i] * vars[i]) + rhs
    CPCONSTRAINT(bool_linear_ge, ~CPVAR(r), CPVAR(ZEROVAR), cp_neg_coeffs, cp_vars, rhs);
}

// Clause Constraints
// ------------------

void nutmeg_cp_bool_clause(SolverInstanceBase& s, const Call* call)
{
    const auto& pos = BOOLVARARRAY(0);
    const auto& neg = BOOLVARARRAY(1);
    vec<geas::clause_elt> clause;
    for (auto& var : pos)
    {
        clause.push(CPVAR(var));
    }
    for (auto& var : neg)
    {
        clause.push(~CPVAR(var));
    }
    CPCLAUSE(clause);
}

// Coercion Constraints
// --------------------

void nutmeg_cp_bool2int(SolverInstanceBase& s, const Call* call)
{
    CPCONSTRAINT(add_clause, CPVAR(BOOLVAR(0)), CPVAR(INTVAR(1)) <= 0);
    CPCONSTRAINT(add_clause, ~CPVAR(BOOLVAR(0)), CPVAR(INTVAR(1)) >= 1);
}
void nutmeg_mip_bool2int(SolverInstanceBase& s, const Call* call) // SPECIAL
{
    SOL.add_mip_int_var_as_bool_var_alias(BOOLVAR(0), INTVAR(1));
}

// Element Constraints
// --------------------

void nutmeg_cp_array_int_element(SolverInstanceBase& s, const Call* call)
{
    release_assert(ARRAY(1)->min(0) == 1 && ARRAY(1)->max(0) == ARRAY(1)->size(),
                   "Error in cp_array_int_element constraint");

    const auto& index_var = INTVAR(0);
    const auto& array = INTARRAY(1);
    const auto& value_var = INTVAR(2);

    vec<int> cp_array;
    for (const auto x : array)
    {
        cp_array.push(x);
    }
    CPCONSTRAINT(int_element, CPVAR(value_var), CPVAR(index_var), cp_array);
}

void nutmeg_cp_array_var_int_element(SolverInstanceBase& s, const Call* call)
{
    release_assert(ARRAY(1)->min(0) == 1 && ARRAY(1)->max(0) == ARRAY(1)->size(),
                   "Error in cp_array_var_int_element constraint");

    const auto& index_var = INTVAR(0);
    const auto& array = INTVARARRAY(1);
    const auto& value_var = INTVAR(2);

    vec<geas::intvar> cp_array;
    for (const auto x : array)
    {
        cp_array.push(CPVAR(x));
    }
    CPCONSTRAINT(var_int_element, CPVAR(value_var), CPVAR(index_var), cp_array);
}

// Global constraints
// ------------------

void nutmeg_cp_all_different(SolverInstanceBase& s, const Call* call)
{
    const auto& vars = INTVARARRAY(0);

    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(all_different_int, cp_vars);
}

void nutmeg_cp_all_different_except_0(SolverInstanceBase& s, const Call* call) {
    const auto& vars = INTVARARRAY(0);

    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_vars.push(CPVAR(vars[i]));
    }
    CPCONSTRAINT(all_different_except_0, cp_vars);
}

void nutmeg_cp_disjunctive(SolverInstanceBase& s, const Call* call)
{
    const auto& start_times = INTVARARRAY(0);
    const auto& durations = INTARRAY(1);

    vec<geas::intvar> cp_start_times;
    vec<int> cp_durations;
    for (size_t i = 0; i < start_times.size(); ++i)
    {
        cp_start_times.push(CPVAR(start_times[i]));
        cp_durations.push(durations[i]);
    }
    CPCONSTRAINT(disjunctive_int, cp_start_times, cp_durations);
}

void nutmeg_cp_cumulative(SolverInstanceBase& s, const Call* call)
{
    const auto& start_times = INTVARARRAY(0);
    const auto& durations = INTARRAY(1);
    const auto& resources = INTARRAY(2);
    const auto capacity = INT(3);

    vec<geas::intvar> cp_start_times;
    vec<int> cp_durations;
    vec<int> cp_resources;
    for (size_t i = 0; i < start_times.size(); ++i)
    {
        cp_start_times.push(CPVAR(start_times[i]));
        cp_durations.push(durations[i]);
        cp_resources.push(resources[i]);
    }
    CPCONSTRAINT(cumulative, cp_start_times, cp_durations, cp_resources, capacity);
}

void nutmeg_cp_cumulative_var(SolverInstanceBase& s, const Call* call)
{
    const auto& start_times = INTVARARRAY(0);
    const auto& durations = INTVARARRAY(1);
    const auto& resources = INTVARARRAY(2);
    const auto capacity = INTVAR(3);

    vec<geas::intvar> cp_start_times;
    vec<geas::intvar> cp_durations;
    vec<geas::intvar> cp_resources;
    for (size_t i = 0; i < start_times.size(); ++i)
    {
        cp_start_times.push(CPVAR(start_times[i]));
        cp_durations.push(CPVAR(durations[i]));
        cp_resources.push(CPVAR(resources[i]));
    }
    CPCONSTRAINT(cumulative_var, cp_start_times, cp_durations, cp_resources, CPVAR(capacity));
}

void nutmeg_cp_cumulative_optional(SolverInstanceBase& s, const Call* call)
{
    const auto& active = BOOLVARARRAY(0);
    const auto& start_times = INTVARARRAY(1);
    const auto& durations = INTVARARRAY(2);
    const auto& resources = INTARRAY(3);
    const auto capacity = INT(4);

    vec<geas::patom_t> cp_active;
    vec<geas::intvar> cp_start_times;
    vec<geas::intvar> cp_durations;
    vec<int> cp_resources;
    for (size_t i = 0; i < start_times.size(); ++i)
    {
        cp_active.push(CPVAR(active[i]));
        cp_start_times.push(CPVAR(start_times[i]));
        cp_durations.push(CPVAR(durations[i]));
        cp_resources.push(resources[i]);
    }
    CPCONSTRAINT(cumulative_sel, cp_start_times, cp_durations, cp_resources, cp_active, capacity);
}

void nutmeg_cp_global_cardinality(SolverInstanceBase& s, const Call* call)
{
    const auto& x = INTVARARRAY(0);
    const auto& cover = INTARRAY(1);
    const auto& count = INTARRAY(2);

    vec<geas::intvar> cp_x;
    vec<int> cp_cover;
    vec<int> cp_count;
    for (size_t i = 0; i < x.size(); ++i)
    {
        cp_x.push(CPVAR(x[i]));
    }
    for (size_t i = 0; i < cover.size(); ++i)
    {
        cp_cover.push(cover[i]);
    }
    for (size_t i = 0; i < count.size(); ++i)
    {
        cp_count.push(count[i]);
    }

    vec<int> srcs(cp_x.size(), 1);
    vec<geas::bflow> flows;
    for (int i = 0; i < cp_x.size(); ++i)
        for (int j = 0; j < cp_cover.size(); ++j)
            if (cp_x[i].lb(CPDATA) <= cp_cover[j] && cp_cover[j] <= cp_x[i].ub(CPDATA))
            {
                flows.push({i, j, cp_x[i] == cp_cover[j]});
            }
    geas::bipartite_flow(CPDATA, srcs, cp_count, flows);
}

void nutmeg_cp_table_int(SolverInstanceBase& s, const Call* call)
{
    const auto& vars = INTVARARRAY(0);
    const auto& tmp = INTARRAY(1);
    release_assert(tmp.size() % vars.size() == 0, "Internal error while adding table constraint to Nutmeg");

    vec<geas::intvar> cp_vars;
    for (size_t i = 0; i < vars.size(); ++i)
    {
        cp_vars.push(CPVAR(vars[i]));
    }

    vec<vec<int>> table(tmp.size() == 0 ? 0 : tmp.size()/vars.size());
    for (int i = 0; i < table.size(); ++i) {
        table[i].growTo(vars.size());
        for (int j = 0; j < vars.size(); ++j) {
            table[i][j] = tmp[i*vars.size() + j];
        }
    }
    geas::table_id id = geas::table::build(CPDATA, table);
    geas::table::post(CPDATA, id, cp_vars);
}

    // void p_bool_or_imp(SolverInstanceBase& s, const Call* call) {
    //   geas::add_clause(SD, ~BOOLVAR(2), BOOLVAR(0), BOOLVAR(1));
    // }

    // void p_bool_and_imp(SolverInstanceBase& s, const Call* call) {
    //   geas::add_clause(SD, ~BOOLVAR(2), BOOLVAR(0));
    //   geas::add_clause(SD, ~BOOLVAR(2), BOOLVAR(1));
    // }

    // void p_bool_xor_imp(SolverInstanceBase& s, const Call* call) {
    //   p_bool_ne_imp(s, call);
    // }

    // void p_bool_clause_imp(SolverInstanceBase& s, const Call* call) {
    //   auto pos = ARRAY(0);
    //   auto neg = ARRAY(1);
    //   vec<geas::clause_elt> clause;
    //   clause.push(~BOOLVAR(2));
    //   for (int i = 0; i < pos->size(); ++i) {
    //     clause.push(SI.asBoolVar((*pos)[i]));
    //   }
    //   for (int j = 0; j < neg->size(); ++j) {
    //     clause.push(~SI.asBoolVar((*neg)[j]));
    //   }
    //   geas::add_clause(*SD, clause);
    // }

    // void p_array_bool_or_imp(SolverInstanceBase& s, const Call* call) {
    //   auto arr = ARRAY(0);
    //   vec<geas::clause_elt> clause;
    //   clause.push(~BOOLVAR(1));
    //   for (int i = 0; i < arr->size(); ++i) {
    //     geas::patom_t elem = SI.asBoolVar((*arr)[i]);
    //     clause.push(elem);
    //   }
    //   geas::add_clause(*SD, clause);
    // }

    // void p_array_bool_and_imp(SolverInstanceBase& s, const Call* call) {
    //   auto arr = ARRAY(0);
    //   for (int i = 0; i < arr->size(); ++i) {
    //     geas::add_clause(SD, ~BOOLVAR(1), SI.asBoolVar((*arr)[i]));
    //   }
    // }

    // void p_bool_clause_reif(SolverInstanceBase& s, const Call* call) {
    //   auto pos = ARRAY(0);
    //   auto neg = ARRAY(1);
    //   vec<geas::clause_elt> clause;
    //   clause.push(~BOOLVAR(2));
    //   for (int i = 0; i < pos->size(); ++i) {
    //     geas::patom_t elem = SI.asBoolVar((*pos)[i]);
    //     geas::add_clause(SD, BOOLVAR(2), ~elem);
    //     clause.push(elem);
    //   }
    //   for (int j = 0; j < neg->size(); ++j) {
    //     geas::patom_t elem = SI.asBoolVar((*neg)[j]);
    //     geas::add_clause(SD, BOOLVAR(2), elem);
    //     clause.push(~elem);
    //   }
    //   geas::add_clause(*SD, clause);
    // }

    // void p_array_bool_element(SolverInstanceBase& s, const Call* call) {
    //   assert(ARRAY(1)->min(0) == 1 && ARRAY(1)->max(0) == ARRAY(1)->size()+1);
    //   vec<bool> vals = BOOLARRAY(1);
    //   if (PAR(0)) {
    //     SOL.post(vals[INT(0)-1] ? BOOLVAR(2) : ~BOOLVAR(2));
    //   } else if (PAR(2)) {
    //     for (int j = 0; j < vals.size(); ++j) {
    //       if (vals[j] != BOOL(2)) {
    //         SOL.post(INTVAR(0) != j+1);
    //       }
    //     }
    //   } else {
    //     for (int j = 0; j < vals.size(); ++j) {
    //       geas::add_clause(SD, INTVAR(0) != j+1, vals[j] ? BOOLVAR(2) : ~BOOLVAR(2));
    //     }
    //   }
    // }

    // void p_array_var_bool_element(SolverInstanceBase& s, const Call* call) {
    //   assert(ARRAY(1)->min(0) == 1 && ARRAY(1)->max(0) == ARRAY(1)->size()+1);
    //   if (PAR(1)) {
    //     return p_array_bool_element(s, call);
    //   }
    //   if (PAR(0) && PAR(2)) {
    //     SOL.post(BOOL(2) ? SI.asBoolVar((*ARRAY(1))[INT(0) - 1]) : ~SI.asBoolVar((*ARRAY(1))[INT(0) - 1]));
    //   } else if (PAR(0)) {
    //     Expression* elem = (*ARRAY(1))[INT(0)-1];
    //     if (elem->type().ispar()) {
    //       return p_array_bool_element(s, call);
    //     } else {
    //       geas::add_clause(SD, BOOLVAR(2), ~SI.asBoolVar(elem));
    //       geas::add_clause(SD, ~BOOLVAR(2), SI.asBoolVar(elem));
    //     }
    //   } else if (PAR(2)) {
    //     for (int j = 0; j < ARRAY(1)->size(); ++j) {
    //       Expression* elem = (*ARRAY(1))[j];
    //       if (elem->type().isvar()) {
    //         geas::add_clause(SD, INTVAR(0) != j+1, INT(2) ? SI.asBoolVar(elem) : ~SI.asBoolVar(elem));
    //       } else {
    //         if (SI.asBool(elem) != INT(2)) {
    //           SOL.post(INTVAR(0) != j+1);
    //         }
    //       }
    //     }
    //   } else {
    //     auto vars = BOOLVARARRAY(1);
    //     for (int j = 0; j < vars.size(); ++j) {
    //       geas::add_clause(SD, INTVAR(0) != j+1, ~vars[j], BOOLVAR(2));
    //       geas::add_clause(SD, INTVAR(0) != j+1, vars[j], ~BOOLVAR(2));
    //     }
    //   }
    // }

    // void p_at_most(SolverInstanceBase& s, const Call* call) {
    //   vec<geas::intvar> ivars = INTVARARRAY(1);
    //   vec<geas::patom_t> bvars;
    //   for (auto &ivar : ivars) {
    //     bvars.push(ivar == INT(2));
    //   }

    //   if (INT(0) == 1) {
    //     geas::atmost_1(SD, bvars);
    //   } else {
    //     geas::atmost_k(SD, bvars, INT(0));
    //   }
    // }

    // void p_at_most1(SolverInstanceBase& s, const Call* call) {
    //   vec<geas::intvar> ivars = INTVARARRAY(0);
    //   vec<geas::patom_t> bvars;
    //   for (auto &ivar : ivars) {
    //     bvars.push(ivar == INT(1));
    //   }
    //   geas::atmost_1(SD, bvars);
    // }

    // void p_disjunctive(SolverInstanceBase& s, const Call* call) {
    //   vec<geas::intvar> st = INTVARARRAY(0);
    //   if (PAR(1)) {
    //     vec<int> d = INTARRAY(1);
    //     geas::disjunctive_int(SD, st, d);
    //   } else {
    //     vec<geas::intvar> d = INTVARARRAY(1);
    //     geas::disjunctive_var(SD, st, d);
    //   }
    // }
  }
}

#pragma clang diagnostic pop