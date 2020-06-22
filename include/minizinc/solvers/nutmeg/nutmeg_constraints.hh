/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Edward Lam <edward.lam@monash.edu>
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <minizinc/ast.hh>
#include <minizinc/solver_instance_base.hh>

namespace MiniZinc {
  namespace NutmegConstraints {

#define PosterImpl(X) void X(SolverInstanceBase& s, const Call* ce)

    /* Integer and Boolean Linear Constraints */
    PosterImpl(nutmeg_cp_int_bool_lin_le);
    PosterImpl(nutmeg_mip_int_bool_lin_le);
    PosterImpl(nutmeg_cp_int_bool_lin_eq);
    PosterImpl(nutmeg_mip_int_bool_lin_eq);

    /* Integer Linear Constraints */
    PosterImpl(nutmeg_cp_int_lin_le);
    PosterImpl(nutmeg_mip_int_lin_le);
    PosterImpl(nutmeg_cp_int_lin_le_imp);
    PosterImpl(nutmeg_mip_int_lin_le_imp);
    PosterImpl(nutmeg_cp_int_lin_le_reif);
    PosterImpl(nutmeg_mip_int_lin_le_reif);
    PosterImpl(nutmeg_cp_int_lin_eq);
    PosterImpl(nutmeg_mip_int_lin_eq);
    PosterImpl(nutmeg_cp_int_lin_eq_imp);
    PosterImpl(nutmeg_mip_int_lin_eq_imp);
    PosterImpl(nutmeg_cp_int_lin_eq_reif);
    PosterImpl(nutmeg_cp_int_lin_ne);
    PosterImpl(nutmeg_cp_int_lin_ne_imp);
    PosterImpl(nutmeg_cp_int_lin_ne_reif);

    /* Integer Comparison Constraints */
    PosterImpl(nutmeg_cp_int_le);
    PosterImpl(nutmeg_cp_int_le_imp);
    PosterImpl(nutmeg_cp_int_le_reif);
    PosterImpl(int_eq);
    PosterImpl(nutmeg_cp_int_eq);
    PosterImpl(nutmeg_mip_int_eq);
    PosterImpl(nutmeg_cp_int_eq_imp);
    PosterImpl(int_eq_reif);
    PosterImpl(nutmeg_cp_int_eq_reif);
    PosterImpl(nutmeg_cp_int_ne);
    PosterImpl(nutmeg_cp_int_ne_imp);
    PosterImpl(nutmeg_cp_int_ne_reif);

    /* Integer Arithmetic Constraints */
    PosterImpl(nutmeg_cp_int_abs);
    PosterImpl(nutmeg_cp_int_times);
    PosterImpl(nutmeg_cp_int_div);
    PosterImpl(nutmeg_cp_int_min);
    PosterImpl(nutmeg_cp_int_max);

    /* Boolean Linear Constraints */
    PosterImpl(nutmeg_cp_bool_lin_le);
    PosterImpl(nutmeg_mip_bool_lin_le);
    PosterImpl(nutmeg_cp_bool_lin_le_imp);
    PosterImpl(nutmeg_mip_bool_lin_le_imp);
    PosterImpl(nutmeg_cp_bool_lin_le_reif);
    PosterImpl(nutmeg_mip_bool_lin_le_reif);
    PosterImpl(nutmeg_cp_bool_lin_eq);
    PosterImpl(nutmeg_mip_bool_lin_eq);
    PosterImpl(nutmeg_cp_bool_lin_eq_imp);
    PosterImpl(nutmeg_mip_bool_lin_eq_imp);
    PosterImpl(nutmeg_cp_bool_lin_eq_reif);
    PosterImpl(nutmeg_cp_bool_lin_ne);
    PosterImpl(nutmeg_cp_bool_lin_ne_imp);
    PosterImpl(nutmeg_cp_bool_lin_ne_reif);

    /* Boolean Arithmetic Constraints */
//    PosterImpl(nutmeg_cp_bool_or_imp);
//    PosterImpl(nutmeg_cp_bool_and_imp);
//    PosterImpl(nutmeg_cp_bool_xor_imp);

    /* Clause Constraints */
    PosterImpl(nutmeg_cp_bool_clause);
//    PosterImpl(nutmeg_cp_bool_clause_imp);
//    PosterImpl(nutmeg_cp_array_bool_or_imp);
//    PosterImpl(nutmeg_cp_array_bool_and_imp);
//    PosterImpl(nutmeg_cp_bool_clause_reif);

    /* Coercion Constraints */
    PosterImpl(nutmeg_cp_bool2int);
    PosterImpl(nutmeg_mip_bool2int);

    /* Element Constraints */
    PosterImpl(nutmeg_cp_array_int_element);
    PosterImpl(nutmeg_cp_array_var_int_element);
//    PosterImpl(nutmeg_cp_array_bool_element);
//    PosterImpl(nutmeg_cp_array_var_bool_element);
//    PosterImpl(nutmeg_cp_array_int_minimum);
//    PosterImpl(nutmeg_cp_array_int_maximum);

    /* Global Constraints */
    PosterImpl(nutmeg_cp_all_different);
    PosterImpl(nutmeg_cp_all_different_except_0);
//    PosterImpl(nutmeg_cp_at_most);
//    PosterImpl(nutmeg_cp_at_most1);
    PosterImpl(nutmeg_cp_disjunctive);
    PosterImpl(nutmeg_cp_cumulative);
    PosterImpl(nutmeg_cp_cumulative_var);
    PosterImpl(nutmeg_cp_cumulative_optional);
    PosterImpl(nutmeg_cp_global_cardinality);
    PosterImpl(nutmeg_cp_table_int);
  }
}
