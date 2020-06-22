include "cumulative.mzn";

% TODO: Propagator for Disjunctive is buggy in Geas

predicate fzn_disjunctive_strict(array[int] of var int: s, array[int] of var int: d) =
    cp_disjunctive_strict(s,d) /\ mip_disjunctive_strict(s,d);

predicate cp_disjunctive_strict(array[int] of var int: s, array[int] of var int: d) =
    if is_fixed(d) then
        nutmeg_cp_cumulative(s, fix(d), [1 | i in index_set(s)], 1)
    else
        nutmeg_cp_cumulative_var(s, d, [1 | i in index_set(s)], 1)
    endif;

predicate mip_disjunctive_strict(array[int] of var int: s, array[int] of var int: d) =
    cp_disjunctive_strict(s,d);