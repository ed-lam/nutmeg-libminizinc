include "nosets.mzn";

% Given a constraint XXX, predicate XXX(...) posts the constraint to both the MIP and CP solvers.
% cp_XXX(...) posts the CP constraint only. mip_XXX(...) posts the MIP constraint (or relaxation) only.
% nutmeg_cp_XXX(...) and nutmeg_mip_XXX(...) are private.

% ----------------------------------------------------------------------------------------------------------------------
% Integer and Boolean Linear Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate int_bool_lin_le(array [int] of int: as, array [int] of var int: bs,
                          array [int] of int: cs, array [int] of var bool: ds,
                          int: e) =
    cp_int_bool_lin_le(as,bs,cs,ds,e) /\ mip_int_bool_lin_le(as,bs,cs,ds,e);
predicate cp_int_bool_lin_le(array [int] of int: as, array [int] of var int: bs,
                             array [int] of int: cs, array [int] of var bool: ds,
                             int: e) =
    if length(as) == 1 then
        let {int: first = min(index_set(as))} in
        if as[first] == 1 then
            nutmeg_cp_int_bool_lin_le(bs[first],cs,ds,e)
        else
            if as[first] == -1 then
                nutmeg_cp_int_bool_lin_le(bs[first],[-cs[i] | i in index_set(cs)],ds,-e)
            else
                let {int: x_min = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then lb(bs[i]) else ub(bs[i]) endif)),
                     int: x_max = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then ub(bs[i]) else lb(bs[i]) endif)),
                     var x_min..x_max: x} in
                cp_int_lin_le(as ++ [-1],bs ++ [x],0) /\ nutmeg_cp_int_bool_lin_le(x,cs,ds,e)
            endif
        endif
    else
        let {int: x_min = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then lb(bs[i]) else ub(bs[i]) endif)),
             int: x_max = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then ub(bs[i]) else lb(bs[i]) endif)),
             var x_min..x_max: x} in
        cp_int_lin_le(as ++ [-1],bs ++ [x],0) /\ nutmeg_cp_int_bool_lin_le(x,cs,ds,e)
    endif;
predicate mip_int_bool_lin_le(array [int] of int: as, array [int] of var int: bs,
                              array [int] of int: cs, array [int] of var bool: ds,
                              int: e) =
    nutmeg_mip_int_bool_lin_le(as,bs,cs,ds,e);
predicate nutmeg_cp_int_bool_lin_le(var int: x,
                                    array [int] of int: cs, array [int] of var bool: ds,
                                    int: e);
predicate nutmeg_mip_int_bool_lin_le(array [int] of int: as, array [int] of var int: bs,
                                     array [int] of int: cs, array [int] of var bool: ds,
                                     int: e);

predicate int_bool_lin_eq(array [int] of int: as, array [int] of var int: bs,
                          array [int] of int: cs, array [int] of var bool: ds,
                          int: e) =
    cp_int_bool_lin_eq(as,bs,cs,ds,e) /\ mip_int_bool_lin_eq(as,bs,cs,ds,e);
predicate cp_int_bool_lin_eq(array [int] of int: as, array [int] of var int: bs,
                             array [int] of int: cs, array [int] of var bool: ds,
                             int: e) =
    if length(as) == 1 then
        let {int: first = min(index_set(as))} in
        if as[first] == 1 then
            nutmeg_cp_int_bool_lin_eq(bs[first],cs,ds,e)
        else
            if as[first] == -1 then
                nutmeg_cp_int_bool_lin_eq(bs[first],[-cs[i] | i in index_set(cs)],ds,-e)
            else
                let {int: x_min = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then lb(bs[i]) else ub(bs[i]) endif)),
                     int: x_max = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then ub(bs[i]) else lb(bs[i]) endif)),
                     var x_min..x_max: x} in
                cp_int_lin_eq(as ++ [-1],bs ++ [x],0) /\ nutmeg_cp_int_bool_lin_eq(x,cs,ds,e)
            endif
        endif
    else
        let {int: x_min = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then lb(bs[i]) else ub(bs[i]) endif)),
             int: x_max = sum(i in index_set(as)) (as[i] * (if as[i] >= 0 then ub(bs[i]) else lb(bs[i]) endif)),
             var x_min..x_max: x} in
        cp_int_lin_eq(as ++ [-1],bs ++ [x],0) /\ nutmeg_cp_int_bool_lin_eq(x,cs,ds,e)
    endif;
predicate mip_int_bool_lin_eq(array [int] of int: as, array [int] of var int: bs,
                              array [int] of int: cs, array [int] of var bool: ds,
                              int: e) =
    nutmeg_mip_int_bool_lin_eq(as,bs,cs,ds,e);
predicate nutmeg_cp_int_bool_lin_eq(var int: x,
                                    array [int] of int: cs, array [int] of var bool: ds,
                                    int: e);
predicate nutmeg_mip_int_bool_lin_eq(array [int] of int: as, array [int] of var int: bs,
                                     array [int] of int: cs, array [int] of var bool: ds,
                                     int: e);

% ----------------------------------------------------------------------------------------------------------------------
% Integer Linear Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate int_lin_lt(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_lt(as,bs,c) /\ mip_int_lin_lt(as,bs,c);
predicate cp_int_lin_lt(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_le(as,bs,c-1);
predicate mip_int_lin_lt(array [int] of int: as, array [int] of var int: bs, int: c) =
    mip_int_lin_le(as,bs,c-1);

predicate int_lin_gt(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_gt(as,bs,c) /\ mip_int_lin_gt(as,bs,c);
predicate cp_int_lin_gt(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_le([-as[i] | i in index_set(as)], bs, -(c+1));
predicate mip_int_lin_gt(array [int] of int: as, array [int] of var int: bs, int: c) =
    mip_int_lin_le([-as[i] | i in index_set(as)], bs, -(c+1));

predicate int_lin_ge(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_ge(as,bs,c) /\ mip_int_lin_ge(as,bs,c);
predicate cp_int_lin_ge(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_le([-as[i] | i in index_set(as)], bs, -c);
predicate mip_int_lin_ge(array [int] of int: as, array [int] of var int: bs, int: c) =
    mip_int_lin_le([-as[i] | i in index_set(as)], bs, -c);

% ----------------------------------------------------------------------------------------------------------------------

predicate int_lin_le(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_le(as,bs,c) /\ mip_int_lin_le(as,bs,c);
predicate cp_int_lin_le(array [int] of int: as, array [int] of var int: bs, int: c) =
    nutmeg_cp_int_lin_le(as,bs,c);
predicate mip_int_lin_le(array [int] of int: as, array [int] of var int: bs, int: c) =
    nutmeg_mip_int_lin_le(as,bs,c);
predicate nutmeg_cp_int_lin_le(array [int] of int: as, array [int] of var int: bs, int: c);
predicate nutmeg_mip_int_lin_le(array [int] of int: as, array [int] of var int: bs, int: c);

predicate int_lin_le_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    cp_int_lin_le_imp(as,bs,c,r) /\ mip_int_lin_le_imp(as,bs,c,r);
predicate cp_int_lin_le_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lin_le(as,bs,c)
        else
            true
        endif
    else
        nutmeg_cp_int_lin_le_imp(as,bs,c,r)
    endif;
predicate mip_int_lin_le_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lin_le(as,bs,c)
        else
            true
        endif
    else
        nutmeg_mip_int_lin_le_imp(as,bs,c,r)
    endif;
predicate nutmeg_cp_int_lin_le_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);
predicate nutmeg_mip_int_lin_le_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);

predicate int_lin_le_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    cp_int_lin_le_reif(as,bs,c,r) /\ mip_int_lin_le_reif(as,bs,c,r);
predicate cp_int_lin_le_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lin_le(as,bs,c)
        else
            cp_int_lin_gt(as,bs,c)
        endif
    else
        nutmeg_cp_int_lin_le_reif(as,bs,c,r)
    endif;
predicate mip_int_lin_le_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lin_le(as,bs,c)
        else
            mip_int_lin_gt(as,bs,c)
        endif
    else
        nutmeg_mip_int_lin_le_reif(as,bs,c,r)
    endif;
predicate nutmeg_cp_int_lin_le_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);
predicate nutmeg_mip_int_lin_le_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);

% ----------------------------------------------------------------------------------------------------------------------

predicate int_lin_eq(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_eq(as,bs,c) /\ mip_int_lin_eq(as,bs,c);
predicate cp_int_lin_eq(array [int] of int: as, array [int] of var int: bs, int: c) =
    nutmeg_cp_int_lin_eq(as,bs,c);
predicate mip_int_lin_eq(array [int] of int: as, array [int] of var int: bs, int: c) =
    nutmeg_mip_int_lin_eq(as,bs,c);
predicate nutmeg_cp_int_lin_eq(array [int] of int: as, array [int] of var int: bs, int: c);
predicate nutmeg_mip_int_lin_eq(array [int] of int: as, array [int] of var int: bs, int: c);

predicate int_lin_eq_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    cp_int_lin_eq_imp(as,bs,c,r) /\ mip_int_lin_eq_imp(as,bs,c,r);
predicate cp_int_lin_eq_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lin_eq(as,bs,c)
        else
            true
        endif
    else
        nutmeg_cp_int_lin_eq_imp(as,bs,c,r)
    endif;
predicate mip_int_lin_eq_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lin_eq(as,bs,c)
        else
            true
        endif
    else
        nutmeg_mip_int_lin_eq_imp(as,bs,c,r)
    endif;
predicate nutmeg_cp_int_lin_eq_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);
predicate nutmeg_mip_int_lin_eq_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);

predicate int_lin_eq_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    cp_int_lin_eq_reif(as,bs,c,r) /\ mip_int_lin_eq_reif(as,bs,c,r);
predicate cp_int_lin_eq_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lin_eq(as,bs,c)
        else
            cp_int_lin_ne(as,bs,c)
        endif
    else
        nutmeg_cp_int_lin_eq_reif(as,bs,c,r)
    endif;
predicate mip_int_lin_eq_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lin_eq(as,bs,c)
        else
            mip_int_lin_ne(as,bs,c)
        endif
    else
        nutmeg_mip_int_lin_eq_reif(as,bs,c,r)
    endif;
predicate nutmeg_cp_int_lin_eq_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);
predicate nutmeg_mip_int_lin_eq_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    mip_int_lin_eq_imp(as,bs,c,r); % RELAXATION

% ----------------------------------------------------------------------------------------------------------------------

predicate int_lin_ne(array [int] of int: as, array [int] of var int: bs, int: c) =
    cp_int_lin_ne(as,bs,c) /\ mip_int_lin_ne(as,bs,c);
predicate cp_int_lin_ne(array [int] of int: as, array [int] of var int: bs, int: c) =
    nutmeg_cp_int_lin_ne(as,bs,c);
predicate mip_int_lin_ne(array [int] of int: as, array [int] of var int: bs, int: c) =
    nutmeg_mip_int_lin_ne(as,bs,c);
predicate nutmeg_cp_int_lin_ne(array [int] of int: as, array [int] of var int: bs, int: c);
predicate nutmeg_mip_int_lin_ne(array [int] of int: as, array [int] of var int: bs, int: c) =
    true; % TODO

predicate int_lin_ne_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    cp_int_lin_ne_imp(as,bs,c,r) /\ mip_int_lin_ne_imp(as,bs,c,r);
predicate cp_int_lin_ne_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lin_ne(as,bs,c)
        else
            true
        endif
    else
        nutmeg_cp_int_lin_ne_imp(as,bs,c,r)
    endif;
predicate mip_int_lin_ne_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lin_ne(as,bs,c)
        else
            true
        endif
    else
        nutmeg_mip_int_lin_ne_imp(as,bs,c,r)
    endif;
predicate nutmeg_cp_int_lin_ne_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);
predicate nutmeg_mip_int_lin_ne_imp(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    true; % TODO

predicate int_lin_ne_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    cp_int_lin_ne_reif(as,bs,c,r) /\ mip_int_lin_ne_reif(as,bs,c,r);
predicate cp_int_lin_ne_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lin_ne(as,bs,c)
        else
            cp_int_lin_eq(as,bs,c)
        endif
    else
        nutmeg_cp_int_lin_ne_reif(as,bs,c,r)
    endif;
predicate mip_int_lin_ne_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lin_ne(as,bs,c)
        else
            mip_int_lin_eq(as,bs,c)
        endif
    else
        nutmeg_mip_int_lin_ne_reif(as,bs,c,r)
    endif;
predicate nutmeg_cp_int_lin_ne_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r);
predicate nutmeg_mip_int_lin_ne_reif(array [int] of int: as, array [int] of var int: bs, int: c, var bool: r) =
    mip_int_lin_eq_imp(as,bs,c,not r); % RELAXATION

% ----------------------------------------------------------------------------------------------------------------------
% Integer Comparison Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate int_gt(var int: a, var int: b) =
    cp_int_gt(a,b) /\ mip_int_gt(a,b);
predicate cp_int_gt(var int: a, var int: b) =
    nutmeg_cp_int_le(b,a,-1);
predicate mip_int_gt(var int: a, var int: b) =
    nutmeg_mip_int_le(b,a,-1);

predicate int_ge(var int: a, var int: b) =
    cp_int_ge(a,b) /\ mip_int_ge(a,b);
predicate cp_int_ge(var int: a, var int: b) =
    nutmeg_cp_int_le(b,a,0);
predicate mip_int_ge(var int: a, var int: b) =
    nutmeg_mip_int_le(b,a,0);

% ----------------------------------------------------------------------------------------------------------------------

predicate int_le(var int: a, var int: b) =
    cp_int_le(a,b) /\ mip_int_le(a,b);
predicate int_le(var int: a, var int: b, int: c) =
    cp_int_le(a,b,c) /\ mip_int_le(a,b,c);
predicate cp_int_le(var int: a, var int: b) =
    cp_int_le(a,b,0);
predicate mip_int_le(var int: a, var int: b) =
    mip_int_le(a,b,0);
predicate cp_int_le(var int: a, var int: b, int: c) =
    nutmeg_cp_int_le(a,b,c);
predicate mip_int_le(var int: a, var int: b, int: c) =
    nutmeg_mip_int_le(a,b,c);
predicate nutmeg_cp_int_le(var int: a, var int: b, int: c);
predicate nutmeg_mip_int_le(var int: a, var int: b, int: c) =
    mip_int_lin_le([1,-1], [a,b], c);

predicate int_le_imp(var int: a, var int: b, var bool: r) =
    cp_int_le_imp(a,b,r) /\ mip_int_le_imp(a,b,r);
predicate int_le_imp(var int: a, var int: b, int: c, var bool: r) =
    cp_int_le_imp(a,b,c,r) /\ mip_int_le_imp(a,b,c,r);
predicate cp_int_le_imp(var int: a, var int: b, var bool: r) =
    cp_int_le_imp(a,b,0,r);
predicate mip_int_le_imp(var int: a, var int: b, var bool: r) =
    mip_int_le_imp(a,b,0,r);
predicate cp_int_le_imp(var int: a, var int: b, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_le(a,b,c)
        else
            true
        endif
    else
        nutmeg_cp_int_le_imp(a,b,c,r)
    endif;
predicate mip_int_le_imp(var int: a, var int: b, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_le(a,b,c)
        else
            true
        endif
    else
        nutmeg_mip_int_le_imp(a,b,c,r)
    endif;
predicate nutmeg_cp_int_le_imp(var int: a, var int: b, int: c, var bool: r);
predicate nutmeg_mip_int_le_imp(var int: a, var int: b, int: c, var bool: r) =
    mip_int_lin_le_imp([1,-1], [a,b], c, r);

predicate int_le_reif(var int: a, var int: b, var bool: r) =
    cp_int_le_reif(a,b,r) /\ mip_int_le_reif(a,b,r);
predicate cp_int_le_reif(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_le(a,b)
        else
            cp_int_gt(a,b)
        endif
    else
        nutmeg_cp_int_le_reif(a,b,0,r)
    endif;
predicate mip_int_le_reif(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_le(a,b)
        else
            mip_int_gt(a,b)
        endif
    else
        nutmeg_mip_int_le_reif(a,b,0,r)
    endif;
predicate nutmeg_cp_int_le_reif(var int: a, var int: b, int: c, var bool: r);
predicate nutmeg_mip_int_le_reif(var int: a, var int: b, int: c, var bool: r) =
    mip_int_lin_le_reif([1,-1], [a,b], c, r);

% ----------------------------------------------------------------------------------------------------------------------

predicate int_lt(var int: a, var int: b) =
    cp_int_lt(a,b) /\ mip_int_lt(a,b);
predicate cp_int_lt(var int: a, var int: b) =
    nutmeg_cp_int_le(a,b,-1);
predicate mip_int_lt(var int: a, var int: b) =
    nutmeg_mip_int_le(a,b,-1);

predicate int_lt_imp(var int: a, var int: b, var bool: r) =
    cp_int_lt_imp(a,b,r) /\ mip_int_lt_imp(a,b,r);
predicate cp_int_lt_imp(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lt(a,b)
        else
            true
        endif
    else
        nutmeg_cp_int_lt_imp(a,b,r)
    endif;
predicate mip_int_lt_imp(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lt(a,b)
        else
            true
        endif
    else
        nutmeg_mip_int_lt_imp(a,b,r)
    endif;
predicate nutmeg_cp_int_lt_imp(var int: a, var int: b, var bool: r) =
    nutmeg_cp_int_le_imp(a,b,-1,r);
predicate nutmeg_mip_int_lt_imp(var int: a, var int: b, var bool: r) =
    nutmeg_mip_int_le_imp(a,b,-1,r);

predicate int_lt_reif(var int: a, var int: b, var bool: r) =
    cp_int_lt_reif(a,b,r) /\ mip_int_lt_reif(a,b,r);
predicate cp_int_lt_reif(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_lt(a,b)
        else
            cp_int_ge(a,b)
        endif
    else
        nutmeg_cp_int_lt_reif(a,b,r)
    endif;
predicate mip_int_lt_reif(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_lt(a,b)
        else
            mip_int_ge(a,b)
        endif
    else
        nutmeg_mip_int_lt_reif(a,b,r)
    endif;
predicate nutmeg_cp_int_lt_reif(var int: a, var int: b, var bool: r) =
    nutmeg_cp_int_le_reif(a,b,-1,r);
predicate nutmeg_mip_int_lt_reif(var int: a, var int: b, var bool: r) =
    nutmeg_mip_int_le_reif(a,b,-1,r);

% ----------------------------------------------------------------------------------------------------------------------

%predicate int_eq(var int: a, var int: b); % TODO: Seems to be a bit problematic if redefining here.
predicate cp_int_eq(var int: a, var int: b) =
    nutmeg_cp_int_eq(a,b);
predicate mip_int_eq(var int: a, var int: b) =
    nutmeg_mip_int_eq(a,b);
predicate nutmeg_cp_int_eq(var int: a, var int: b);
predicate nutmeg_mip_int_eq(var int: a, var int: b);

predicate int_eq_imp(var int: a, var int: b, var bool: r) =
    cp_int_eq_imp(a,b,r) /\ mip_int_eq_imp(a,b,r);
predicate cp_int_eq_imp(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_eq(a,b)
        else
            true
        endif
    else
        nutmeg_cp_int_eq_imp(a,b,r)
    endif;
predicate mip_int_eq_imp(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_eq(a,b)
        else
            true
        endif
    else
        nutmeg_mip_int_eq_imp(a,b,r)
    endif;
predicate nutmeg_cp_int_eq_imp(var int: a, var int: b, var bool: r);
predicate nutmeg_mip_int_eq_imp(var int: a, var int: b, var bool: r) =
    mip_int_lin_eq_imp([1,-1], [a,b], 0, r);

% No need for int_eq_reif. It is handled by pulling literals out from the solvers.

% ----------------------------------------------------------------------------------------------------------------------

predicate int_ne(var int: a, var int: b) =
    cp_int_ne(a,b) /\ mip_int_ne(a,b);
predicate cp_int_ne(var int: a, var int: b) =
    nutmeg_cp_int_ne(a,b);
predicate mip_int_ne(var int: a, var int: b) =
    nutmeg_mip_int_ne(a,b);
predicate nutmeg_cp_int_ne(var int: a, var int: b);
predicate nutmeg_mip_int_ne(var int: a, var int: b) =
    true; % TODO

predicate int_ne_imp(var int: a, var int: b, var bool: r) =
    cp_int_ne_imp(a,b,r) /\ mip_int_ne_imp(a,b,r);
predicate cp_int_ne_imp(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_ne(a,b)
        else
            true
        endif
    else
        nutmeg_cp_int_ne_imp(a,b,r)
    endif;
predicate mip_int_ne_imp(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_ne(a,b)
        else
            true
        endif
    else
        nutmeg_mip_int_ne_imp(a,b,r)
    endif;
predicate nutmeg_cp_int_ne_imp(var int: a, var int: b, var bool: r);
predicate nutmeg_mip_int_ne_imp(var int: a, var int: b, var bool: r) =
    true; % TODO

predicate int_ne_reif(var int: a, var int: b, var bool: r) =
    cp_int_ne_reif(a,b,r) /\ mip_int_ne_reif(a,b,r);
predicate cp_int_ne_reif(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_int_ne(a,b)
        else
            cp_int_eq(a,b)
        endif
    else
        nutmeg_cp_int_ne_reif(a,b,r)
    endif;
predicate mip_int_ne_reif(var int: a, var int: b, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_int_ne(a,b)
        else
            mip_int_eq(a,b)
        endif
    else
        nutmeg_mip_int_ne_reif(a,b,r)
    endif;
predicate nutmeg_cp_int_ne_reif(var int: a, var int: b, var bool: r);
predicate nutmeg_mip_int_ne_reif(var int: a, var int: b, var bool: r) =
    true; % TODO

% ----------------------------------------------------------------------------------------------------------------------
% Integer Arithmetic Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate int_abs(var int: a, var int: b) =
    cp_int_abs(a,b) /\ mip_int_abs(a,b);
predicate cp_int_abs(var int: a, var int: b) =
    nutmeg_cp_int_abs(a,b);
predicate mip_int_abs(var int: a, var int: b) =
    nutmeg_mip_int_abs(a,b);
predicate nutmeg_cp_int_abs(var int: a, var int: b);
predicate nutmeg_mip_int_abs(var int: a, var int: b) =
    if lb(a) >= 0 then
        mip_int_eq(b,a)
    elseif ub(a) <= 0 then
        mip_int_eq(b,-a)
    else
        mip_int_ge(b,a) /\
        mip_int_ge(b,-a) /\
        mip_int_ge(b,0) /\
        mip_int_le(b, max([ub(a), -lb(a)])) % RELAXATION
    endif;

predicate int_times(var int: a, var int: b, var int: c) =
    cp_int_times(a,b,c) /\ mip_int_times(a,b,c);
predicate cp_int_times(var int: a, var int: b, var int: c) =
    nutmeg_cp_int_times(a,b,c);
predicate mip_int_times(var int: a, var int: b, var int: c) =
    nutmeg_mip_int_times(a,b,c);
predicate nutmeg_cp_int_times(var int: a, var int: b, var int: c);
predicate nutmeg_mip_int_times(var int: a, var int: b, var int: c) =
    true; % TODO

predicate int_div(var int: a, var int: b, var int: c) =
    cp_int_div(a,b,c) /\ mip_int_div(a,b,c);
predicate cp_int_div(var int: a, var int: b, var int: c) =
    nutmeg_cp_int_div(a,b,c);
predicate mip_int_div(var int: a, var int: b, var int: c) =
    nutmeg_mip_int_div(a,b,c);
predicate nutmeg_cp_int_div(var int: a, var int: b, var int: c);
predicate nutmeg_mip_int_div(var int: a, var int: b, var int: c) =
    true; % TODO

predicate int_min(var int: a, var int: b, var int: c) =
    cp_int_min(a,b,c) /\ mip_int_min(a,b,c);
predicate cp_int_min(var int: a, var int: b, var int: c) =
    nutmeg_cp_int_min(a,b,c);
predicate mip_int_min(var int: a, var int: b, var int: c) =
    nutmeg_mip_int_min(a,b,c);
predicate nutmeg_cp_int_min(var int: a, var int: b, var int: c);
predicate nutmeg_mip_int_min(var int: a, var int: b, var int: c) =
    mip_int_le(c,a) /\ mip_int_le(c,b); % RELAXATION

predicate int_max(var int: a, var int: b, var int: c) =
    cp_int_max(a,b,c) /\ mip_int_max(a,b,c);
predicate cp_int_max(var int: a, var int: b, var int: c) =
    nutmeg_cp_int_max(a,b,c);
predicate mip_int_max(var int: a, var int: b, var int: c) =
    nutmeg_mip_int_max(a,b,c);
predicate nutmeg_cp_int_max(var int: a, var int: b, var int: c);
predicate nutmeg_mip_int_max(var int: a, var int: b, var int: c) =
    mip_int_ge(c,a) /\ mip_int_ge(c,b); % RELAXATION

% ----------------------------------------------------------------------------------------------------------------------
% Boolean Linear Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate bool_lin_lt(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_lt(as,bs,c) /\ mip_bool_lin_lt(as,bs,c);
predicate cp_bool_lin_lt(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_le(as,bs,c-1);
predicate mip_bool_lin_lt(array [int] of int: as, array [int] of var bool: bs, int: c) =
    mip_bool_lin_le(as,bs,c-1);

predicate bool_lin_gt(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_gt(as,bs,c) /\ mip_bool_lin_gt(as,bs,c);
predicate cp_bool_lin_gt(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_le([-as[i] | i in index_set(as)], bs, -(c+1));
predicate mip_bool_lin_gt(array [int] of int: as, array [int] of var bool: bs, int: c) =
    mip_bool_lin_le([-as[i] | i in index_set(as)], bs, -(c+1));

predicate bool_lin_ge(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_ge(as,bs,c) /\ mip_bool_lin_ge(as,bs,c);
predicate cp_bool_lin_ge(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_le([-as[i] | i in index_set(as)], bs, -c);
predicate mip_bool_lin_ge(array [int] of int: as, array [int] of var bool: bs, int: c) =
    mip_bool_lin_le([-as[i] | i in index_set(as)], bs, -c);

% ----------------------------------------------------------------------------------------------------------------------

predicate bool_lin_le(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_le(as,bs,c) /\ mip_bool_lin_le(as,bs,c);
predicate cp_bool_lin_le(array [int] of int: as, array [int] of var bool: bs, int: c) =
    nutmeg_cp_bool_lin_le(as,bs,c);
predicate mip_bool_lin_le(array [int] of int: as, array [int] of var bool: bs, int: c) =
    nutmeg_mip_bool_lin_le(as,bs,c);
predicate nutmeg_cp_bool_lin_le(array [int] of int: as, array [int] of var bool: bs, int: c);
predicate nutmeg_mip_bool_lin_le(array [int] of int: as, array [int] of var bool: bs, int: c);

predicate bool_lin_le_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    cp_bool_lin_le_imp(as,bs,c,r) /\ mip_bool_lin_le_imp(as,bs,c,r);
predicate cp_bool_lin_le_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_bool_lin_le(as,bs,c)
        else
            true
        endif
    else
        nutmeg_cp_bool_lin_le_imp(as,bs,c,r)
    endif;
predicate mip_bool_lin_le_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_bool_lin_le(as,bs,c)
        else
            true
        endif
    else
        nutmeg_mip_bool_lin_le_imp(as,bs,c,r)
    endif;
predicate nutmeg_cp_bool_lin_le_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);
predicate nutmeg_mip_bool_lin_le_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);

predicate bool_lin_le_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    cp_bool_lin_le_reif(as,bs,c,r) /\ mip_bool_lin_le_reif(as,bs,c,r);
predicate cp_bool_lin_le_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_bool_lin_le(as,bs,c)
        else
            cp_bool_lin_gt(as,bs,c)
        endif
    else
        nutmeg_cp_bool_lin_le_reif(as,bs,c,r)
    endif;
predicate mip_bool_lin_le_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_bool_lin_le(as,bs,c)
        else
            mip_bool_lin_gt(as,bs,c)
        endif
    else
        nutmeg_mip_bool_lin_le_reif(as,bs,c,r)
    endif;
predicate nutmeg_cp_bool_lin_le_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);
predicate nutmeg_mip_bool_lin_le_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);

% ----------------------------------------------------------------------------------------------------------------------

predicate bool_lin_eq(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_eq(as,bs,c) /\ mip_bool_lin_eq(as,bs,c);
predicate cp_bool_lin_eq(array [int] of int: as, array [int] of var bool: bs, int: c) =
    nutmeg_cp_bool_lin_eq(as,bs,c);
predicate mip_bool_lin_eq(array [int] of int: as, array [int] of var bool: bs, int: c) =
    nutmeg_mip_bool_lin_eq(as,bs,c);
predicate nutmeg_cp_bool_lin_eq(array [int] of int: as, array [int] of var bool: bs, int: c);
predicate nutmeg_mip_bool_lin_eq(array [int] of int: as, array [int] of var bool: bs, int: c);

predicate bool_lin_eq_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    cp_bool_lin_eq_imp(as,bs,c,r) /\ mip_bool_lin_eq_imp(as,bs,c,r);
predicate cp_bool_lin_eq_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_bool_lin_eq(as,bs,c)
        else
            true
        endif
    else
        nutmeg_cp_bool_lin_eq_imp(as,bs,c,r)
    endif;
predicate mip_bool_lin_eq_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_bool_lin_eq(as,bs,c)
        else
            true
        endif
    else
        nutmeg_mip_bool_lin_eq_imp(as,bs,c,r)
    endif;
predicate nutmeg_cp_bool_lin_eq_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);
predicate nutmeg_mip_bool_lin_eq_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);

predicate bool_lin_eq_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    cp_bool_lin_eq_reif(as,bs,c,r) /\ mip_bool_lin_eq_reif(as,bs,c,r);
predicate cp_bool_lin_eq_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            cp_bool_lin_eq(as,bs,c)
        else
            assert(false, "Cannot use bool_lin_eq_reif because bool_lin_ne is not yet implemented in Geas\n") /\
            cp_bool_lin_ne(as,bs,c)
        endif
    else
        assert(false, "Cannot use bool_lin_eq_reif because bool_lin_ne is not yet implemented in Geas\n") /\
        nutmeg_cp_bool_lin_eq_reif(as,bs,c,r)
    endif;
predicate mip_bool_lin_eq_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            mip_bool_lin_eq(as,bs,c)
        else
            assert(false, "Cannot use bool_lin_eq_reif because bool_lin_ne is not yet implemented in Geas\n") /\
            mip_bool_lin_ne(as,bs,c)
        endif
    else
        assert(false, "Cannot use bool_lin_eq_reif because bool_lin_ne is not yet implemented in Geas\n") /\
        nutmeg_mip_bool_lin_eq_reif(as,bs,c,r)
    endif;
predicate nutmeg_cp_bool_lin_eq_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);
predicate nutmeg_mip_bool_lin_eq_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    mip_bool_lin_eq_imp(as,bs,c,r); % RELAXATION

% ----------------------------------------------------------------------------------------------------------------------

predicate bool_lin_ne(array [int] of int: as, array [int] of var bool: bs, int: c) =
    cp_bool_lin_ne(as,bs,c) /\ mip_bool_lin_ne(as,bs,c);
predicate cp_bool_lin_ne(array [int] of int: as, array [int] of var bool: bs, int: c) =
    assert(false, "Cannot use bool_lin_ne because it is not yet implemented in Geas\n") /\
    nutmeg_cp_bool_lin_ne(as,bs,c);
predicate mip_bool_lin_ne(array [int] of int: as, array [int] of var bool: bs, int: c) =
    assert(false, "Cannot use bool_lin_ne because it is not yet implemented in Geas\n") /\
    nutmeg_mip_bool_lin_ne(as,bs,c);
predicate nutmeg_cp_bool_lin_ne(array [int] of int: as, array [int] of var bool: bs, int: c);
predicate nutmeg_mip_bool_lin_ne(array [int] of int: as, array [int] of var bool: bs, int: c) =
    true; % TODO

predicate bool_lin_ne_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    cp_bool_lin_ne_imp(as,bs,c,r) /\ mip_bool_lin_ne_imp(as,bs,c,r);
predicate cp_bool_lin_ne_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            assert(false, "Cannot use bool_lin_ne_imp because bool_lin_ne is not yet implemented in Geas\n") /\
            cp_bool_lin_ne(as,bs,c)
        else
            true
        endif
    else
        assert(false, "Cannot use bool_lin_ne_imp because bool_lin_ne is not yet implemented in Geas\n") /\
        nutmeg_cp_bool_lin_ne_imp(as,bs,c,r)
    endif;
predicate mip_bool_lin_ne_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            assert(false, "Cannot use bool_lin_ne_imp because bool_lin_ne is not yet implemented in Geas\n") /\
            mip_bool_lin_ne(as,bs,c)
        else
            true
        endif
    else
        assert(false, "Cannot use bool_lin_ne_imp because bool_lin_ne is not yet implemented in Geas\n") /\
        nutmeg_mip_bool_lin_ne_imp(as,bs,c,r)
    endif;
predicate nutmeg_cp_bool_lin_ne_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);
predicate nutmeg_mip_bool_lin_ne_imp(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    true; % TODO

predicate bool_lin_ne_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    cp_bool_lin_ne_reif(as,bs,c,r) /\ mip_bool_lin_ne_reif(as,bs,c,r);
predicate cp_bool_lin_ne_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            assert(false, "Cannot use bool_lin_ne_reif because bool_lin_ne is not yet implemented in Geas\n") /\
            cp_bool_lin_ne(as,bs,c)
        else
            cp_bool_lin_eq(as,bs,c)
        endif
    else
        assert(false, "Cannot use bool_lin_ne_reif because bool_lin_ne is not yet implemented in Geas\n") /\
        nutmeg_cp_bool_lin_ne_reif(as,bs,c,r)
    endif;
predicate mip_bool_lin_ne_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            assert(false, "Cannot use bool_lin_ne_reif because bool_lin_ne is not yet implemented in Geas\n") /\
            mip_bool_lin_ne(as,bs,c)
        else
            mip_bool_lin_eq(as,bs,c)
        endif
    else
        assert(false, "Cannot use bool_lin_ne_reif because bool_lin_ne is not yet implemented in Geas\n") /\
        nutmeg_mip_bool_lin_ne_reif(as,bs,c,r)
    endif;
predicate nutmeg_cp_bool_lin_ne_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r);
predicate nutmeg_mip_bool_lin_ne_reif(array [int] of int: as, array [int] of var bool: bs, int: c, var bool: r) =
    mip_bool_lin_eq_imp(as,bs,c,not r); % RELAXATION

% ----------------------------------------------------------------------------------------------------------------------
% Boolean Comparison Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate bool_gt(var bool: a, var bool: b) =
    bool_lt(b,a);

predicate bool_ge(var bool: a, var bool: b) =
    cp_bool_ge(a,b) /\ mip_bool_ge(a,b);
predicate cp_bool_ge(var bool: a, var bool: b) =
    cp_bool_le(b,a);
predicate mip_bool_ge(var bool: a, var bool: b) =
    mip_bool_le(b,a);

% ----------------------------------------------------------------------------------------------------------------------

% a <= b
% ~b -> ~a
% b \/ ~a
predicate bool_le(var bool: a, var bool: b) =
    cp_bool_le(a,b) /\ mip_bool_le(a,b);
predicate cp_bool_le(var bool: a, var bool: b) =
    nutmeg_cp_bool_le(a,b);
predicate mip_bool_le(var bool: a, var bool: b) =
    nutmeg_mip_bool_le(a,b);
predicate nutmeg_cp_bool_le(var bool: a, var bool: b) =
    cp_bool_clause([b],[a]);
predicate nutmeg_mip_bool_le(var bool: a, var bool: b) =
    mip_bool_clause([b],[a]);

% r -> a <= b
% r -> (~b -> ~a)
% r -> b \/ ~a
% ~r \/ b \/ ~a
predicate bool_le_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_le_imp(a,b,r) /\ mip_bool_le_imp(a,b,r);
predicate cp_bool_le_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_le_imp(a,b,r);
predicate mip_bool_le_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_le_imp(a,b,r);
predicate nutmeg_cp_bool_le_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_clause([b],[a,r]);
predicate nutmeg_mip_bool_le_imp(var bool: a, var bool: b, var bool: r) =
    mip_bool_clause([b],[a,r]);

% Above r -> a <= b and:
% ~r -> a > b
% ~r -> a /\ ~b
% (~r -> a) /\ (~r -> ~b)
% (r \/ a) /\ (r \/ ~b)
predicate bool_le_reif(var bool: a, var bool: b, var bool: r) =
    cp_bool_le_reif(a,b,r) /\ mip_bool_le_reif(a,b,r);
predicate cp_bool_le_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_le_reif(a,b,r);
predicate mip_bool_le_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_le_reif(a,b,r);
predicate nutmeg_cp_bool_le_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_le_imp(a,b,r) /\ cp_bool_clause([r,a],[]) /\ cp_bool_clause([r],[b]);
predicate nutmeg_mip_bool_le_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_le_imp(a,b,r) /\ mip_bool_clause([r,a],[]) /\ mip_bool_clause([r],[b]);

% ----------------------------------------------------------------------------------------------------------------------

% a < b
% ~a /\ b
predicate bool_lt(var bool: a, var bool: b) =
    (not a) /\ b;

% r -> a < b
% r -> ~a /\ b
% ~r \/ (~a /\ b)
% (~r \/ ~a) /\ (~r \/ b)
predicate bool_lt_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_lt_imp(a,b,r) /\ mip_bool_lt_imp(a,b,r);
predicate cp_bool_lt_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_lt_imp(a,b,r);
predicate mip_bool_lt_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_lt_imp(a,b,r);
predicate nutmeg_cp_bool_lt_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_clause([],[r,a]) /\ cp_bool_clause([b],[r]);
predicate nutmeg_mip_bool_lt_imp(var bool: a, var bool: b, var bool: r) =
    mip_bool_clause([],[r,a]) /\ mip_bool_clause([b],[r]);

% Above r -> a < b and:
% ~r -> a >= b
% r \/ (~a -> ~b)
% r \/ (a \/ ~b)
% r \/ a \/ ~b
predicate bool_lt_reif(var bool: a, var bool: b, var bool: r) =
    cp_bool_lt_reif(a,b,r) /\ mip_bool_lt_reif(a,b,r);
predicate cp_bool_lt_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_lt_reif(a,b,r);
predicate mip_bool_lt_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_lt_reif(a,b,r);
predicate nutmeg_cp_bool_lt_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_lt_imp(a,b,r) /\ cp_bool_clause([r,a],[b]);
predicate nutmeg_mip_bool_lt_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_lt_imp(a,b,r) /\ mip_bool_clause([r,a],[b]);

% ----------------------------------------------------------------------------------------------------------------------

% a == b
% (a -> b) /\ (~a -> ~b)
% (~a \/ b) /\ (a \/ ~b)
predicate bool_eq(var bool: a, var bool: b) =
    cp_bool_eq(a,b) /\ mip_bool_eq(a,b);
predicate cp_bool_eq(var bool: a, var bool: b) =
    nutmeg_cp_bool_eq(a,b);
predicate mip_bool_eq(var bool: a, var bool: b) =
    nutmeg_mip_bool_eq(a,b);
predicate nutmeg_cp_bool_eq(var bool: a, var bool: b) =
    cp_bool_clause([b],[a]) /\ cp_bool_clause([a],[b]);
predicate nutmeg_mip_bool_eq(var bool: a, var bool: b) =
    mip_bool_clause([b],[a]) /\ mip_bool_clause([a],[b]);

% r -> a == b
% (r -> (a -> b)) /\ (r -> (~a -> ~b))
% (~r \/ (~a \/ b)) /\ (~r \/ (a \/ ~b))
% (~r \/ ~a \/ b) /\ (~r \/ a \/ ~b)
predicate bool_eq_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_eq_imp(a,b,r) /\ mip_bool_eq_imp(a,b,r);
predicate cp_bool_eq_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_eq_imp(a,b,r);
predicate mip_bool_eq_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_eq_imp(a,b,r);
predicate nutmeg_cp_bool_eq_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_clause([b],[r,a]) /\ cp_bool_clause([a],[r,b]);
predicate nutmeg_mip_bool_eq_imp(var bool: a, var bool: b, var bool: r) =
    mip_bool_clause([b],[r,a]) /\ mip_bool_clause([a],[r,b]);

% Above r -> a == b and:
% ~r -> a != b
% (r \/ ~a \/ ~b) /\ (r \/ a \/ b) - copied from bool_ne_imp
predicate bool_eq_reif(var bool: a, var bool: b, var bool: r) =
    cp_bool_eq_reif(a,b,r) /\ mip_bool_eq_reif(a,b,r);
predicate cp_bool_eq_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_eq_reif(a,b,r);
predicate mip_bool_eq_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_eq_reif(a,b,r);
predicate nutmeg_cp_bool_eq_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_eq_imp(a,b,r) /\ cp_bool_clause([r],[a,b]) /\ cp_bool_clause([r,a,b],[]);
predicate nutmeg_mip_bool_eq_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_eq_imp(a,b,r) /\ mip_bool_clause([r],[a,b]) /\ mip_bool_clause([r,a,b],[]);

% ----------------------------------------------------------------------------------------------------------------------

% a != b
% (a -> ~b) /\ (~a -> b)
% (~a \/ ~b) /\ (a \/ b)
predicate bool_ne(var bool: a, var bool: b) =
    cp_bool_ne(a,b) /\ mip_bool_ne(a,b);
predicate cp_bool_ne(var bool: a, var bool: b) =
    nutmeg_cp_bool_ne(a,b);
predicate mip_bool_ne(var bool: a, var bool: b) =
    nutmeg_mip_bool_ne(a,b);
predicate nutmeg_cp_bool_ne(var bool: a, var bool: b) =
    cp_bool_clause([],[a,b]) /\ cp_bool_clause([a,b],[]);
predicate nutmeg_mip_bool_ne(var bool: a, var bool: b) =
    mip_bool_clause([],[a,b]) /\ mip_bool_clause([a,b],[]);

% r -> a != b
% r -> (a -> ~b) /\ (~a -> b)
% ~r \/ (~a \/ ~b) /\ (a \/ b)
% (~r \/ ~a \/ ~b) /\ (~r \/ a \/ b)
predicate bool_ne_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_ne_imp(a,b,r) /\ mip_bool_ne_imp(a,b,r);
predicate cp_bool_ne_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_ne_imp(a,b,r);
predicate mip_bool_ne_imp(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_ne_imp(a,b,r);
predicate nutmeg_cp_bool_ne_imp(var bool: a, var bool: b, var bool: r) =
    cp_bool_clause([],[r,a,b]) /\ cp_bool_clause([a,b],[r]);
predicate nutmeg_mip_bool_ne_imp(var bool: a, var bool: b, var bool: r) =
    mip_bool_clause([],[r,a,b]) /\ mip_bool_clause([a,b],[r]);

% Above r -> a != b and:
% ~r -> a == b
% (r \/ ~a \/ b) /\ (r \/ a \/ ~b) - copied from bool_eq_imp
predicate bool_ne_reif(var bool: a, var bool: b, var bool: r) =
    cp_bool_ne_reif(a,b,r) /\ mip_bool_ne_reif(a,b,r);
predicate cp_bool_ne_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_ne_reif(a,b,r);
predicate mip_bool_ne_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_ne_reif(a,b,r);
predicate nutmeg_cp_bool_ne_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_ne_imp(a,b,r) /\ cp_bool_clause([r,b],[a]) /\ cp_bool_clause([r,a],[b]);
predicate nutmeg_mip_bool_ne_reif(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_ne_imp(a,b,r) /\ mip_bool_clause([r,b],[a]) /\ mip_bool_clause([r,a],[b]);

% ----------------------------------------------------------------------------------------------------------------------
% Boolean Arithmetic Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate bool_or(var bool: a, var bool: b, var bool: r) =
    cp_bool_or(a,b,r) /\ mip_bool_or(a,b,r);
predicate cp_bool_or(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_or(a,b,r);
predicate mip_bool_or(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_or(a,b,r);
predicate nutmeg_cp_bool_or(var bool: a, var bool: b, var bool: r) =
    cp_array_bool_or([a,b],r);
predicate nutmeg_mip_bool_or(var bool: a, var bool: b, var bool: r) =
    mip_array_bool_or([a,b],r);

predicate bool_and(var bool: a, var bool: b, var bool: r) =
    cp_bool_and(a,b,r) /\ mip_bool_and(a,b,r);
predicate cp_bool_and(var bool: a, var bool: b, var bool: r) =
    nutmeg_cp_bool_and(a,b,r);
predicate mip_bool_and(var bool: a, var bool: b, var bool: r) =
    nutmeg_mip_bool_and(a,b,r);
predicate nutmeg_cp_bool_and(var bool: a, var bool: b, var bool: r) =
    cp_array_bool_and([a,b],r);
predicate nutmeg_mip_bool_and(var bool: a, var bool: b, var bool: r) =
    mip_array_bool_and([a,b],r);

predicate bool_xor(var bool: a, var bool: b) =
    bool_ne(a,b);

predicate bool_xor(var bool: a, var bool: b, var bool: r) =
    bool_ne_reif(a,b,r);

predicate bool_not(var bool: a, var bool: b) =
    bool_ne(a,b);

% ----------------------------------------------------------------------------------------------------------------------
% Clause Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate bool_clause(array [int] of var bool: as,array [int] of var bool: bs) =
    cp_bool_clause(as,bs) /\ mip_bool_clause(as,bs);
predicate cp_bool_clause(array [int] of var bool: as,array [int] of var bool: bs) =
    nutmeg_cp_bool_clause(as,bs);
predicate mip_bool_clause(array [int] of var bool: as,array [int] of var bool: bs) =
    nutmeg_mip_bool_clause(as,bs);
predicate nutmeg_cp_bool_clause(array [int] of var bool: as,array [int] of var bool: bs);
predicate nutmeg_mip_bool_clause(array [int] of var bool: as,array [int] of var bool: bs) =
    mip_bool_lin_le([-1 | i in index_set(as)] ++ [1 | i in index_set(bs)],
                    as ++ bs,
                    -1 + length(bs));

predicate array_bool_or(array [int] of var bool: as, var bool: r) =
    cp_array_bool_or(as,r) /\ mip_array_bool_or(as,r);
predicate cp_array_bool_or(array [int] of var bool: as, var bool: r) =
    nutmeg_cp_array_bool_or(as,r);
predicate mip_array_bool_or(array [int] of var bool: as, var bool: r) =
    nutmeg_mip_array_bool_or(as,r);
predicate nutmeg_cp_array_bool_or(array [int] of var bool: as, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            % r -> a1 \/ ... \/ an
            % ~r \/ (a1 \/ ... \/ an)
            % a1 \/ ... \/ an \/ ~r
            cp_bool_clause(as,[r])
        else
            % ~r -> ~(a1 \/ ... \/ an)
            % r \/ ~(a1 \/ ... \/ an)
            % r \/ (~a1 /\ ... /\ ~an)
            % (r \/ ~a1) /\ ... /\ (r \/ ~an)
            forall(i in index_set(as)) (cp_bool_clause([r],[as[i]]))
        endif
    else
        cp_bool_clause(as,[r]) /\ forall(i in index_set(as)) (cp_bool_clause([r],[as[i]]))
    endif;
predicate nutmeg_mip_array_bool_or(array [int] of var bool: as, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            % r -> a1 \/ ... \/ an
            % ~r \/ (a1 \/ ... \/ an)
            % a1 \/ ... \/ an \/ ~r
            mip_bool_clause(as,[r])
        else
            % ~r -> ~(a1 \/ ... \/ an)
            % r \/ ~(a1 \/ ... \/ an)
            % r \/ (~a1 /\ ... /\ ~an)
            % (r \/ ~a1) /\ ... /\ (r \/ ~an)
            forall(i in index_set(as)) (mip_bool_clause([r],[as[i]]))
        endif
    else
        mip_bool_clause(as,[r]) /\ forall(i in index_set(as)) (mip_bool_clause([r],[as[i]]))
    endif;

predicate array_bool_and(array [int] of var bool: as, var bool: r) =
    cp_array_bool_and(as,r) /\ mip_array_bool_and(as,r);
predicate cp_array_bool_and(array [int] of var bool: as, var bool: r) =
    nutmeg_cp_array_bool_and(as,r);
predicate mip_array_bool_and(array [int] of var bool: as, var bool: r) =
    nutmeg_mip_array_bool_and(as,r);
predicate nutmeg_cp_array_bool_and(array [int] of var bool: as, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            % r -> a1 /\ ... /\ an
            % ~r \/ (a1 /\ ... /\ an)
            % (~r \/ a1) /\ ... /\ (~r \/ an)
            % (a1 \/ ~r) /\ ... /\ (an \/ ~r)
            forall(i in index_set(as)) (cp_bool_clause([as[i]],[r]))
        else
            % ~r -> ~(a1 /\ ... /\ an)
            % r \/ ~(a1 /\ ... /\ an)
            % r \/ ~a1 \/ ... \/ ~an
            cp_bool_clause([r],as)
        endif
    else
        forall(i in index_set(as)) (cp_bool_clause([as[i]],[r])) /\ cp_bool_clause([r],as)
    endif;
predicate nutmeg_mip_array_bool_and(array [int] of var bool: as, var bool: r) =
    if is_fixed(r) then
        if fix(r) then
            % r -> a1 /\ ... /\ an
            % ~r \/ (a1 /\ ... /\ an)
            % (~r \/ a1) /\ ... /\ (~r \/ an)
            % (a1 \/ ~r) /\ ... /\ (an \/ ~r)
            forall(i in index_set(as)) (mip_bool_clause([as[i]],[r]))
        else
            % ~r -> ~(a1 /\ ... /\ an)
            % r \/ ~(a1 /\ ... /\ an)
            % r \/ ~a1 \/ ... \/ ~an
            mip_bool_clause([r],as)
        endif
    else
        forall(i in index_set(as)) (mip_bool_clause([as[i]],[r])) /\ mip_bool_clause([r],as)
    endif;

% ----------------------------------------------------------------------------------------------------------------------
% Coercion Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate bool2int(var bool: a, var int: b) =
    cp_bool2int(a,b) /\ mip_bool2int(a,b);
predicate cp_bool2int(var bool: a, var int: b) =
    nutmeg_cp_bool2int(a,b);
predicate mip_bool2int(var bool: a, var int: b) =
    nutmeg_mip_bool2int(a,b);
predicate nutmeg_cp_bool2int(var bool: a, var int: b);
predicate nutmeg_mip_bool2int(var bool: a, var int: b);

% ----------------------------------------------------------------------------------------------------------------------
% Element Constraints
% ----------------------------------------------------------------------------------------------------------------------

% b == i -> c == as[i]
predicate array_int_element(var int: b, array [int] of int: as, var int: c) =
    cp_array_int_element(b,as,c) /\ mip_array_int_element(b,as,c);
predicate cp_array_int_element(var int: b, array [int] of int: as, var int: c) =
    b in index_set(as) /\ nutmeg_cp_array_int_element(b,as,c);
predicate mip_array_int_element(var int: b, array [int] of int: as, var int: c) =
    b in index_set(as) /\ nutmeg_mip_array_int_element(b,as,c);
predicate nutmeg_cp_array_int_element(var int: b, array [int] of int: as, var int: c);
predicate nutmeg_mip_array_int_element(var int: b, array [int] of int: as, var int: c) =
    forall(i in index_set(as)) (mip_int_lin_eq_imp([1], [c], as[i], (b == i)));

% b == i -> c == as[i]
predicate array_var_int_element(var int: b, array [int] of var int: as, var int: c) =
    cp_array_var_int_element(b,as,c) /\ mip_array_var_int_element(b,as,c);
predicate cp_array_var_int_element(var int: b, array [int] of var int: as, var int: c) =
    b in index_set(as) /\ nutmeg_cp_array_var_int_element(b,as,c);
predicate mip_array_var_int_element(var int: b, array [int] of var int: as, var int: c) =
    b in index_set(as) /\ nutmeg_mip_array_var_int_element(b,as,c);
predicate nutmeg_cp_array_var_int_element(var int: b, array [int] of var int: as, var int: c);
predicate nutmeg_mip_array_var_int_element(var int: b, array [int] of var int: as, var int: c) =
    forall(i in index_set(as)) (mip_int_eq_imp(c, as[i], (b == i)));

% b == i -> c == as[i]
predicate array_bool_element(var int: b, array [int] of bool: as, var bool: c) =
    cp_array_bool_element(b,as,c) /\ mip_array_bool_element(b,as,c);
predicate cp_array_bool_element(var int: b, array [int] of bool: as, var bool: c) =
    b in index_set(as) /\ nutmeg_cp_array_bool_element(b,as,c);
predicate mip_array_bool_element(var int: b, array [int] of bool: as, var bool: c) =
    b in index_set(as) /\ nutmeg_mip_array_bool_element(b,as,c);
predicate nutmeg_cp_array_bool_element(var int: b, array [int] of bool: as, var bool: c) =
    assert(min(index_set(as)) == 1, "Array in array_bool_element must be indexed starting from 1") /\
    % (b == i) -> c == as[i]
    % ~(b == i) \/ c == as[i]
    forall(i in index_set(as)) (cp_bool_clause([if as[i] then c else not c endif], [b == i]));
predicate nutmeg_mip_array_bool_element(var int: b, array [int] of bool: as, var bool: c) =
    forall(i in index_set(as)) (mip_bool_lin_eq_imp([1], [c], as[i], (b == i)));

% b == i -> c == as[i]
predicate array_var_bool_element(var int: b, array [int] of var bool: as, var bool: c) =
    cp_array_var_bool_element(b,as,c) /\ mip_array_var_bool_element(b,as,c);
predicate cp_array_var_bool_element(var int: b, array [int] of var bool: as, var bool: c) =
    b in index_set(as) /\ nutmeg_cp_array_var_bool_element(b,as,c);
predicate mip_array_var_bool_element(var int: b, array [int] of var bool: as, var bool: c) =
    b in index_set(as) /\ nutmeg_mip_array_var_bool_element(b,as,c);
predicate nutmeg_cp_array_var_bool_element(var int: b, array [int] of var bool: as, var bool: c) =
    assert(min(index_set(as)) == 1, "Array in array_var_bool_element must be indexed starting from 1") /\
    % (b == i) -> (as[i] -> c) /\ (~as[i] -> ~c)
    % ~(b == i) \/ ((~as[i] \/ c) /\ (as[i] \/ ~c))
    % (~(b == i) \/ ~as[i] \/ c) /\ (~(b == i) \/ as[i] \/ ~c)
    forall(i in index_set(as)) (cp_bool_clause([c], [b == i, as[i]])) /\
    forall(i in index_set(as)) (cp_bool_clause([as[i]], [b == i, c]));
predicate nutmeg_mip_array_var_bool_element(var int: b, array [int] of var bool: as, var bool: c) =
    forall(i in index_set(as)) (mip_bool_eq_imp(c, as[i], (b == i)));

% ----------------------------------------------------------------------------------------------------------------------
% Set Constraints
% ----------------------------------------------------------------------------------------------------------------------

predicate set_in(var int: x, set of int: s) =
    cp_set_in(x,s) /\ mip_set_in(x,s);
predicate cp_set_in(var int: x, set of int: s) =
    nutmeg_cp_set_in(x,s);
predicate mip_set_in(var int: x, set of int: s) =
    nutmeg_mip_set_in(x,s);
predicate nutmeg_cp_set_in(var int: x, set of int: s) =
    cp_bool_clause([x == k | k in s],[]);
predicate nutmeg_mip_set_in(var int: x, set of int: s) =
    forall(k in lb(x)..ub(x) where not (k in s)) (mip_bool_clause([], [x == k])); % RELAXATION

predicate set_in_reif(var int: x, set of int: s, var bool: r) =
    cp_set_in_reif(x,s,r) /\ mip_set_in_reif(x,s,r);
predicate cp_set_in_reif(var int: x, set of int: s, var bool: r) =
    nutmeg_cp_set_in_reif(x,s,r);
predicate mip_set_in_reif(var int: x, set of int: s, var bool: r) =
    nutmeg_mip_set_in_reif(x,s,r);
predicate nutmeg_cp_set_in_reif(var int: x, set of int: s, var bool: r) =
    %if card(s) = max(s) - min(s) + 1 then
    %    r <-> (x >= min(s) /\ x <= max(s))
    %else
        cp_bool_clause([x == k | k in s], [r]) /\
        forall(k in s) (cp_bool_clause([r], [x == k]));
    %endif; % TODO
predicate nutmeg_mip_set_in_reif(var int: x, set of int: s, var bool: r) =
    % r -> min(s) <= x <= max(s)
    mip_int_lin_le_imp([1], [x], max(s), r) /\
    mip_int_lin_le_imp([-1], [x], -min(s), r); % RELAXATION
