function define_getATC!(mod::Model)
    TCONNECT = mod.ext[:parameters][:TCONNECT]
    n_int    = length(TCONNECT)


    # all sign combinations for the redispatch poly-tope
    signs = collect(Iterators.product(fill((-1, 1), n_int)...))

    mod.ext[:sets][:JV] = 1:length(signs)     # vertices of the box
    mod.ext[:parameters][:signs] = signs
    JN = mod.ext[:sets][:JN]

    # Get node symbols from common parameters
    node_syms = mod.ext[:parameters][:nodes]

    node_shares = Dict{String, Dict{Symbol,Float64}}()
    for zone in zones
        cons_conf = data["Consumers"][zone]
        raw_ns = get(cons_conf, "NodeShare", Dict())
        node_shares[zone] = Dict(Symbol(k) => v for (k,v) in raw_ns)
    end
    
    # Node share matrix: zones × nodes
    M = [get(node_shares[zone], node, 0.0) for zone in zones, node in node_syms]

    mod.ext[:parameters][:CapCM_nodal] = vec(M' * mod.ext[:parameters][:CapCM_zonal])
    return mod
end