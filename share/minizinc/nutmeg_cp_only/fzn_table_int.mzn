predicate fzn_table_int(array [int] of var int: x, array [int, int] of int: t) =
    cp_table_int(x, [ t[i,j] | i in index_set_1of2(t), j in index_set_2of2(t) ])
    /\
    mip_table_int(x, [ t[i,j] | i in index_set_1of2(t), j in index_set_2of2(t) ]);

predicate cp_table_int(array [int] of var int: x, array [int] of int: t) =
    nutmeg_cp_table_int(x,t);
predicate mip_table_int(array [int] of var int: x, array [int] of int: t) =
    cp_table_int(x,t);

predicate nutmeg_cp_table_int(array [int] of var int: x, array [int] of int: t);
