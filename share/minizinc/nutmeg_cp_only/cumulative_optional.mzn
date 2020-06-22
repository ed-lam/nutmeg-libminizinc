include "fzn_cumulative_optional.mzn";

predicate cumulative_optional(array[int] of var bool: active,
                              array[int] of var int: start,
                              array[int] of int: duration,
                              array[int] of int: resource,
                              int: capacity) =
    assert(index_set(active) == index_set(start) /\
           index_set(active) == index_set(duration) /\
           index_set(active) == index_set(resource),
           "cumulative_optional: the 4 array arguments must have identical index sets") /\
    assert(lb_array(duration) >= 0 /\ lb_array(resource) >= 0,
           "cumulative_optional: durations and resource usages must be non-negative") /\
    fzn_cumulative_optional(active, start, duration, resource, capacity);

predicate cumulative_optional(array[int] of var bool: active,
                              array[int] of var int: start,
                              array[int] of int: duration,
                              array[int] of int: resource,
                              int: capacity,
                              var int: makespan) =
    assert(index_set(active) == index_set(start) /\
           index_set(active) == index_set(duration) /\
           index_set(active) == index_set(resource),
           "cumulative_optional: the 4 array arguments must have identical index sets") /\
    assert(lb_array(duration) >= 0 /\ lb_array(resource) >= 0,
           "cumulative_optional: durations and resource usages must be non-negative") /\
    fzn_cumulative_optional(active, start, duration, resource, capacity, makespan);
