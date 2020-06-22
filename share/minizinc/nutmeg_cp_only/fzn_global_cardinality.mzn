predicate fzn_global_cardinality(array[int] of var int: x, array[int] of int: cover, array[int] of var int: count) =
    cp_global_cardinality(x, cover, count) /\ mip_global_cardinality(x, cover, count);

predicate cp_global_cardinality(array[int] of var int: x, array[int] of int: cover, array[int] of var int: count) =
    if is_fixed(count) then
        nutmeg_cp_global_cardinality(x,cover,fix(count))
    else
        assert(false, "global_cardinality constraint for variable count array is not yet implemented")
    endif;
predicate mip_global_cardinality(array[int] of var int: x, array[int] of int: cover, array[int] of var int: count) =
    cp_global_cardinality(x,cover,count);

predicate nutmeg_cp_global_cardinality(array[int] of var int: x, array[int] of int: cover, array[int] of int: count);

%predicate fzn_global_cardinality(array[int] of var int: x, array[int] of int: cover, array[int] of var int: count) =
%    forall (i in index_set(cover)) (
%        count[i] = sum (j in index_set(x)) (bool2int(x[j] = cover[i]))
%    );