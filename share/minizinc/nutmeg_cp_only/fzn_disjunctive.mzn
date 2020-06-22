include "cumulative.mzn";

% TODO: Propagator for Disjunctive is buggy in Geas

predicate fzn_disjunctive(array[int] of var int: s, array[int] of var int: d) =
    cp_disjunctive(s,d) /\ mip_disjunctive(s,d);

predicate cp_disjunctive(array[int] of var int: s, array[int] of var int: d) =
    if is_fixed(d) then
        nutmeg_cp_cumulative(s, fix(d), [1 | i in index_set(s)], 1)
    else
        nutmeg_cp_cumulative_var(s, d, [1 | i in index_set(s)], 1)
    endif;
predicate mip_disjunctive(array[int] of var int: s, array[int] of var int: d) =
    cp_disjunctive(s,d);