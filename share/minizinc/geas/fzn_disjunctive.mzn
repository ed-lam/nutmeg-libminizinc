%predicate geas_disjunctive_var(array[int] of var int: s,
%  array[int] of var int: d);

%predicate disjunctive(array[int] of var int: s, array[int] of var int: d) =
%  geas_disjunctive_var(s, d);
predicate fzn_disjunctive(array[int] of var int: s, array[int] of int: d) =
  geas_cumulative(s, d, [1 | i in index_set(s)], 1);
