predicate fzn_cumulative(array[int] of var int: s, array[int] of int: d, array[int] of int: r, int: b) =
    cp_cumulative(s, d, r, b) /\ mip_cumulative(s, d, r, b);

predicate cp_cumulative(array[int] of var int: s, array[int] of int: d, array[int] of int: r, int: b) =
    nutmeg_cp_cumulative(s,d,r,b);
predicate mip_cumulative(array[int] of var int: s, array[int] of int: d, array[int] of int: r, int: b) =
    cp_cumulative(s,d,r,b);

predicate nutmeg_cp_cumulative(array[int] of var int: s, array[int] of int: d, array[int] of int: r, int: b);

predicate fzn_cumulative(array[int] of var int: s, array[int] of var int: d, array[int] of var int: r, var int: b) =
    cp_cumulative_var(s, d, r, b) /\ mip_cumulative_var(s, d, r, b);

predicate cp_cumulative_var(array[int] of var int: s, array[int] of var int: d, array[int] of var int: r, var int: b) =
    nutmeg_cp_cumulative_var(s,d,r,b);
predicate mip_cumulative_var(array[int] of var int: s, array[int] of var int: d, array[int] of var int: r, var int: b) =
    cp_cumulative_var(s,d,r,b);

predicate nutmeg_cp_cumulative_var(array[int] of var int: s, array[int] of var int: d, array[int] of var int: r, var int: b);