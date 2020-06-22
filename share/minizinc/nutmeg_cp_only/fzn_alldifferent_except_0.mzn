predicate fzn_alldifferent_except_0(array [int] of var int: x) =
    cp_all_different_except_0(x) /\ mip_all_different_except_0(x);

predicate cp_all_different_except_0(array [int] of var int: x) =
    nutmeg_cp_all_different_except_0(x);
predicate mip_all_different_except_0(array [int] of var int: x) =
    cp_all_different_except_0(x);

predicate nutmeg_cp_all_different_except_0(array [int] of var int: x);