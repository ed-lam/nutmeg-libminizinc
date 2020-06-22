predicate fzn_table_int(array [int] of var int: x, array [int, int] of int: t) =
    cp_table_int(x,t)
    /\
    mip_table_int(x,t);

predicate cp_table_int(array [int] of var int: x, array [int, int] of int: t) =
    nutmeg_cp_table_int(x, [t[i,j] | i in index_set_1of2(t), j in index_set_2of2(t)]);
predicate mip_table_int(array [int] of var int: x, array [int, int] of int: t) =
    nutmeg_mip_table_int(x,t);

predicate nutmeg_cp_table_int(array [int] of var int: x, array [int] of int: t);
predicate nutmeg_mip_table_int(array [int] of var int: x, array [int, int] of int: t) =
    true;
    % This linearisation is slower than none at all
    %let { set of int: it = index_set_1of2(t), array[it] of var bool: lambda } in
    %trace("table\n") /\
    %mip_bool_lin_eq([1 | i in it], [lambda[i] | i in it], 1)
    %/\
    %forall(j in index_set(x)) (
    %    mip_int_bool_lin_eq([-1], [x[j]],
    %                        [t[i,j] | i in it], [lambda[i] | i in it],
    %                        0)
    %);