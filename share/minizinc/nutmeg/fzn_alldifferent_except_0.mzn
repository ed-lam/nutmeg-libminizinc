predicate fzn_alldifferent_except_0(array [int] of var int: x) =
    cp_all_different_except_0(x) /\ mip_all_different_except_0(x);

predicate cp_all_different_except_0(array [int] of var int: x) =
    nutmeg_cp_all_different_except_0(x);
predicate mip_all_different_except_0(array [int] of var int: x) =
    nutmeg_mip_all_different_except_0(x);

predicate nutmeg_cp_all_different_except_0(array [int] of var int: x);
predicate nutmeg_mip_all_different_except_0(array [int] of var int: x) =
    if length(x) <= 1 then
        true
    else
        let {
            int: min_dom = min([min(dom(x[i])) | i in index_set(x)]),
            int: max_dom = max([max(dom(x[i])) | i in index_set(x)]),
            array[index_set(x),min_dom..max_dom] of var bool: inds = array2d(
                index_set(x),
                min_dom..max_dom,
                [if j in dom(x[i]) then (x[i] == j) else false endif |
                 i in index_set(x), j in min_dom..max_dom]
            );
        } in (
            %trace("  all_different_except_0: x[" ++ show(index_set(x)) ++ "]\n") /\
            forall(d in min_dom..max_dom where d != 0) (
                mip_bool_lin_le([1 | i in index_set(x)], [inds[i,d] | i in index_set(x)], 1)
                % sum(i in index_set(x)) (inds[i,d]) <= 1
            )
        )
    endif;