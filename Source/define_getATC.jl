function define_getATC!(mod::Model)

    TCONNECT = mod.ext[:parameters][:TCONNECT]
    num_INT  = length(TCONNECT)

    # build the 2^num_INT ±1‐sign vectors
    signs = collect(Iterators.product(fill((-1,1), num_INT)...))
    mod.ext[:sets][:JV]         = 1:length(signs)
    mod.ext[:parameters][:signs] = signs

    return mod
end