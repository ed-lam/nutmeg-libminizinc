include "alldifferent_except_0.mzn";

predicate all_different_except_0(array[$X] of var int: vs) =
   alldifferent_except_0(array1d(vs));