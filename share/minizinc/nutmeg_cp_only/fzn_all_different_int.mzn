predicate fzn_all_different_int(array [int] of var int: x) =
    cp_all_different(x) /\ mip_all_different(x);

predicate cp_all_different(array [int] of var int: x) =
    nutmeg_cp_all_different(x);
predicate mip_all_different(array [int] of var int: x) =
    cp_all_different(x);

predicate nutmeg_cp_all_different(array [int] of var int: x);