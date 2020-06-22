predicate fzn_cumulative_optional(array[int] of var bool: active,
                                  array[int] of var int: start,
                                  array[int] of int: duration,
                                  array[int] of int: resource,
                                  int: capacity) =
    nutmeg_cp_cumulative_optional(active, start, duration, resource, capacity);

predicate fzn_cumulative_optional(array[int] of var bool: active,
                                  array[int] of var int: start,
                                  array[int] of int: duration,
                                  array[int] of int: resource,
                                  int: capacity,
                                  var int: makespan) =
    %forall(t in index_set(active)) (active[t] -> makespan >= start[t] + duration[t]) /\ % TODO: Maybe put this in CP only
    forall(t in index_set(active)) (cp_int_le_imp(start[t], makespan, -duration[t], active[t])) /\
    nutmeg_cp_cumulative_optional(active, start, duration, resource, capacity) /\
    nutmeg_mip_cumulative_optional(active, start, duration, resource, capacity, makespan);

predicate nutmeg_cp_cumulative_optional(array[int] of var bool: active,
                                        array[int] of var int: start,
                                        array[int] of int: duration,
                                        array[int] of int: resource,
                                        int: capacity);
predicate nutmeg_mip_cumulative_optional(array[int] of var bool: active,
                                         array[int] of var int: start,
                                         array[int] of int: duration,
                                         array[int] of int: resource,
                                         int: capacity,
                                         var int: makespan) =
    % sum(i in index_set(active)) (resource[i] * duration[i] * active[i]) <= capacity * makespan
    mip_int_bool_lin_le([-capacity], [makespan],
                        [resource[i] * duration[i] | i in index_set(active)], active,
                        0); % RELAXATION;
