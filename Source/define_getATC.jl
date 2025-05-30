function define_getATC!(mod::Model)
    TCONNECT = mod.ext[:parameters][:TCONNECT]
    n_int    = length(TCONNECT)

    # all sign combinations for the redispatch poly-tope
    signs = collect(Iterators.product(fill((-1, 1), n_int)...))

    mod.ext[:sets][:JV] = 1:length(signs)     # vertices of the box
    mod.ext[:parameters][:signs] = signs
    JN = mod.ext[:sets][:JN]
    # mod.ext[:parameters][:getATC_demand] = zeros(length(JN))
    mod.ext[:parameters][:getATC_demand] = fill(0.0, length(JN))
    return mod
end