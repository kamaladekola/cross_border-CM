

# import Pkg
using JuMP
using Gurobi
using DataFrames, CSV, YAML
import MathOptInterface as MOI
using Printf

# Helper functions
function get_any(d::Dict, keys::AbstractVector{<:AbstractString})
    for k in keys
        if haskey(d, k);            return d[k]            end
        if haskey(d, lowercase(k)); return d[lowercase(k)] end
        if haskey(d, uppercase(k)); return d[uppercase(k)] end
    end
    error("Missing keys $(join(keys, ", ")) in dictionary.")
end
get_any(d::Dict, key::AbstractString) = get_any(d, [key])

function safe_col(df::DataFrame, name::Union{String,Symbol})
    s = Symbol(name)
    @assert s ∈ propertynames(df) "Missing column $(s) in $(nameof(parentmodule(df))) DataFrame"
    return df[!, s]
end

struct Sets
H::Int
Z::Int
N::Int
I::Int
L::Int
JH::UnitRange{Int}
JZ::UnitRange{Int}
JN::UnitRange{Int}
JI::UnitRange{Int}
JL::UnitRange{Int}
zones::Vector{String}
nodes::Vector{Symbol}
techs::Vector{Symbol}
end

struct Maps
    zone_of_node::Vector{Int} # length N, zone index per node
    zone_to_idx::Dict{String,Int}
    tech_to_idx::Dict{Symbol,Int}
    zone_nodes_idx::Vector{Vector{Int}} # per z, list of node indices
end

struct Params
    W::Matrix{Float64} # weights per t × z
    D_zonal::Matrix{Float64} # H × Z (demand)
    AV::Array{Float64,3} # H × I × N (availability factor)
    y_node::Matrix{Float64} # I × N (existing capacity)
    MC::Matrix{Float64} # I × Z (marginal cost)
    max_cap::Matrix{Float64} # I × Z (maximum capacity)
    A::Matrix{Float64} # I × Z (availability)
    IC::Matrix{Float64} # I × Z (investment cost)
    S::Matrix{Float64} # Z × N (demand share per zone-node)
    PTDF::Matrix{Float64} # L × N (power transfer distribution factors)
    Fmax::Vector{Float64} # L (thermal limits)
end

function read_inputs(base_dir::AbstractString)
    data = YAML.load_file(joinpath(base_dir, "Input", "config.yaml"))
    load = CSV.read(joinpath(base_dir, "Input", "load.csv"), DataFrame; delim = ";")
    pv = CSV.read(joinpath(base_dir, "Input", "pv.csv"), DataFrame; delim = ";")
    wind_on = CSV.read(joinpath(base_dir, "Input", "wind_onshore.csv"), DataFrame; delim = ";")
    nptdf = CSV.read(joinpath(base_dir, "Input", "nodal_ptdf.csv"), DataFrame; delim = ";")
    lines = CSV.read(joinpath(base_dir, "Input", "lines.csv"), DataFrame; delim = ";")

    H = Int(data["General"]["nTimesteps"])
    zones = collect_zones(data)
    Z = length(zones)
    
    weights = joinpath(base_dir, "Input", "weights.csv")
    W = if isfile(weights)
        df = CSV.read(weights, DataFrame; delim = ";")
        # Read zone-specific columns or use uniform weights
        W_matrix = zeros(Float64, H, Z)
        for (zidx, z) in enumerate(zones)
            col_sym = Symbol(z)  # Changed from Symbol("(z)") to Symbol(z)
            if col_sym ∈ propertynames(df)
                W_matrix[:, zidx] = Float64.(df[!, col_sym][1:H])
            # elseif :weights ∈ propertynames(df)
            #     # Fallback: use single weights column for all zones
            #     W_matrix[:, zidx] = Float64.(df[!, :weights][1:H])
            else
                error("No weight column found for zone $(z). Expected '$(z)' or 'weights' in weights.csv")
            end
        end
        W_matrix
    else
        error("weights.csv file not found")  # Also fixed the error message here
    end

    return (data=data, load=load, pv=pv, wind_on=wind_on,
    nptdf=nptdf, lines=lines, weights=W)
end

function collect_zones(data)
    return sort!(String.(collect(keys(data["Consumers"]))))
end

function collect_techs(data)

    tech_names = String[]
    for (_zone, gens) in data["Generators"]
        for k in keys(gens)
            push!(tech_names, String(k))
        end
    end
    tech_names = unique(tech_names)
    sort!(tech_names)
    return Symbol.(tech_names)
end

function build_sets_and_maps(data, lines::DataFrame)
    zones = collect_zones(data)
    techs = collect_techs(data)
    zone_nodes = Dict{String,Vector{Symbol}}(z => Symbol.(data["Network"]["ZoneMap"][z]) for z in zones)
    nodes = unique(vcat([Symbol.(data["Network"]["ZoneMap"][z]) for z in zones]...))


    H = Int(data["General"]["nTimesteps"]) ; Z = length(zones) ; N = length(nodes)
    I = length(techs) ; L = nrow(lines)


    zone_to_idx = Dict(z => i for (i, z) in enumerate(zones))
    tech_to_idx = Dict(techs[i] => i for i in 1:I)

    zone_of_node = similar(Vector{Int}(), N)
    resize!(zone_of_node, N)

    for (nidx, n) in enumerate(nodes)
        found = false
        for (zi, z) in enumerate(zones)
            if n ∈ zone_nodes[z]
                zone_of_node[nidx] = zi
                found = true
                break
            end
        end
        @assert found "Node $(n) not found in any zone"
    end


    zone_nodes_idx = [Int[] for _ in 1:Z]
    for nidx in 1:N
        z = zone_of_node[nidx]
        push!(zone_nodes_idx[z], nidx)
    end


    sets = Sets(H, Z, N, I, L, 1:H, 1:Z, 1:N, 1:I, 1:L, zones, nodes, techs)
    maps = Maps(zone_of_node, zone_to_idx, tech_to_idx, zone_nodes_idx)
    return sets, maps
end


function define_parameters(data, load, pv, wind_on, nptdf, lines, weights, sets::Sets, maps::Maps)
    H, Z, N, I, L = sets.H, sets.Z, sets.N, sets.I, sets.L

    # demand per zone
    D_zonal = Array{Float64}(undef, H, Z)
    for (zidx, z) in enumerate(sets.zones)
        D_zonal[:, zidx] = Float64.(safe_col(load, z))[1:H]
    end


    # elastic/price parameters per zone
    ela = zeros(Float64, Z)
    WTP = zeros(Float64, Z)
    for (zidx, z) in enumerate(sets.zones)
        ela[zidx] = Float64(get_any(data["Consumers"][z], "ela"))
        WTP[zidx] = Float64(get_any(data["Consumers"][z], "WTP"))
    end

    # demand share per zone-node (Z × N)
    S = zeros(Float64, Z, N)
    for (zidx, z) in enumerate(sets.zones)
        local_nodes = maps.zone_nodes_idx[zidx]
        if haskey(data["Consumers"][z], "NodeShare")
            raw = Dict(Symbol(k)=>Float64(v) for (k,v) in data["Consumers"][z]["NodeShare"])
            tot = sum(get(raw, sets.nodes[n], 0.0) for n in local_nodes)
            if tot == 0
                for n in local_nodes; S[zidx, n] = 1.0/length(local_nodes); end
            else
                for n in local_nodes; S[zidx, n] = get(raw, sets.nodes[n], 0.0) / tot; end
            end
        else
            for n in local_nodes; S[zidx, n] = 1.0/length(local_nodes); end
        end
        @assert isapprox(sum(S[zidx, local_nodes]), 1.0; atol=1e-10) "Demand shares for zone $(z) must sum to 1"
    end

    # availability factors H × I × N
    function zone_af(tech::Symbol, z::String)
        if tech in (:PV, :WindOnshore)
            df = tech === :PV ? pv : wind_on
            return Float64.(safe_col(df, z))[1:H]
        else
            return ones(Float64, H)
        end
    end
    
    AV = zeros(Float64, H, I, N)
    for (i, tech) in enumerate(sets.techs), (zidx, z) in enumerate(sets.zones)
        af = zone_af(tech, z)
        for nidx in maps.zone_nodes_idx[zidx]
            @inbounds AV[:, i, nidx] .= af
        end
    end


    # existing capacity per (i,n)
    y_node = zeros(Float64, I, N)
    for (zidx, z) in enumerate(sets.zones)
        local_nodes = maps.zone_nodes_idx[zidx]
        for (tname, genco) in data["Generators"][z]
            i = maps.tech_to_idx[Symbol(tname)]
            y0 = Float64(get_any(genco, "C"))
            if haskey(genco, "NodeShare")
                shares = Dict(Symbol(k)=>Float64(v) for (k,v) in genco["NodeShare"])
                total = sum(get(shares, sets.nodes[j], 0.0) for j in local_nodes)
                total = total == 0 ? 1.0 : total
                for n in local_nodes
                    s = get(shares, sets.nodes[n], 0.0) / total
                    y_node[i, n] += s * y0
                end
            else
                for n in local_nodes
                    y_node[i, n] += (1.0/length(local_nodes)) * y0
                end
            end
        end
    end


    # variable & investment costs per (i,z)
    MC = zeros(Float64, I, Z)
    A = zeros(Float64, I, Z)
    IC = zeros(Float64, I, Z)
    max_cap = zeros(Float64, I, Z)
    for (zidx, z) in enumerate(sets.zones), (i, tech) in enumerate(sets.techs)
        if haskey(data["Generators"][z], String(tech))
            genco = data["Generators"][z][String(tech)]
        else
            error("Tech $(tech) missing for zone $(z). Add it or provide a default; current code assumes every (i,z) exists.")
        end

        MC[i, zidx] = Float64(get_any(genco, "b"))
        A[i, zidx] = Float64(get_any(genco, "a"))
        IC[i, zidx] = Float64(get_any(genco, "I"))
        max_cap[i, zidx] = Float64(get_any(genco, "max_cap"))
    end

    # PTDF + limits
    PTDF = Matrix{Float64}(undef, L, N)
    for (j, n) in enumerate(sets.nodes)
        PTDF[:, j] = Float64.(safe_col(nptdf, String(n)))
    end
    # println(DataFrame(PTDF, :auto))

    Fmax = Float64.(safe_col(lines, :Fmax))
    @assert size(PTDF,1) == L "PTDF rows must match number of lines"

    return Params(weights, D_zonal, AV, y_node, MC, max_cap, A, IC, S, PTDF, Fmax), ela, WTP
end

function build_planner_eom(; data, load, pv, wind_on, nodal_ptdf_df, lines, weights)
    sets, maps = build_sets_and_maps(data, lines)
    params, ela, WTP = define_parameters(data, load, pv, wind_on, nodal_ptdf_df, lines, weights, sets, maps)

    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag"=>0))

    H, Z, N, I, L = sets.H, sets.Z, sets.N, sets.I, sets.L
    JH, JZ, JN, JI, JL = sets.JH, sets.JZ, sets.JN, sets.JI, sets.JL
    W, D_zonal, AV, y_node = params.W, params.D_zonal, params.AV, params.y_node
    MC, max_cap, A, IC, PTDF, Fmax, S = params.MC, params.max_cap, params.A, params.IC, params.PTDF, params.Fmax, params.S
    zone_of_node = maps.zone_of_node
   
    # Variables
    @variable(model, y[JI, JZ] >= 0)
    @variable(model, y_bar[JI, JN] >= 0)
    @variable(model, g[JH, JI, JZ] >= 0)
    @variable(model, g_bar[JH, JI, JN] >= 0)
    @variable(model, r[JH, JN])
    @variable(model, f[JH, JL])
    @variable(model, p[JH, JZ])
    @variable(model, d_inel[JH, JZ] >= 0)
    @variable(model, d_ela[JH, JZ]  >= 0)
    @variable(model, ens[JH, JZ] >= 0)

    # negative utility function per zone
    function consumer_utility(t, z)
        w = WTP[z]
        e = ela[z]
        if e == 0.0 || D_zonal[t, z] == 0.0
            return -w * d_inel[t, z]
        else
            return -w * (d_inel[t, z] + d_ela[t, z]) + (w / (2 * e * D_zonal[t, z])) * d_ela[t, z]^2
        end
    end


    # marginal cost + quadratic cost
    @expression(model, gen_cost,
        sum(W[t,z] * MC[i, z] * g[t, i, z] for t in JH, i in JI, z in JZ)
        + sum(W[t,z] * A[i, z] / 2 * g[t, i, z]^2 for t in JH, i in JI, z in JZ)
    )
    # investment cost
    @expression(model, inv_cost, sum(IC[i, z] * y[i, z] for i in JI, z in JZ))

    # consumer utility
    @expression(model, neg_utility,  sum(W[t,z] * consumer_utility(t, z) for t in JH, z in JZ))

    # Objective: minimize (costs - utility)
    @objective(model, Min, (inv_cost + gen_cost + neg_utility))

    # Constraints

    # demand
    @constraint(model, [t in JH, z in JZ], d_inel[t, z] == (1 - ela[z]) * D_zonal[t, z] - ens[t, z])
    @constraint(model, [t in JH, z in JZ], d_ela[t, z]  <=  ela[z] * D_zonal[t, z])


    # zonal balance constraints
    bal = @constraint(model, bal[t in JH, z in JZ],  sum(g[t, i, z] for i in JI) - p[t, z] - (d_inel[t, z] + d_ela[t, z]) == 0)

    # Aggregate nodal net positions to zonal net positions
    agg = @constraint(model, [t in JH, z in JZ], 
        p[t, z] == sum(r[t, n] for n in JN if zone_of_node[n] == z))

    # Nodal balance constraints
    nbal = @constraint(model, [t in JH, n in JN], 
        r[t, n] == sum(g_bar[t, i, n] for i in JI) - 
                   (d_inel[t, zone_of_node[n]] + d_ela[t, zone_of_node[n]]) * S[zone_of_node[n], n])

    # nodal capacity limits
    cap  = @constraint(model, [t in JH, i in JI, n in JN], g_bar[t, i, n] <= AV[t, i, n] * (y_node[i, n] + y_bar[i, n] ))

    # Mainly constrains renewables investment
    Ren_cap = @constraint(model, [i in JI, z in JZ], sum(y_node[i, n] for n in JN if zone_of_node[n] == z) + y[i, z] <= max_cap[i, z])

    # link between zonal and nodal capacity (distorts zonal prices)
    alloc = @constraint(model, [i in JI, z in JZ], y[i, z] == sum(y_bar[i, n] for n in JN if zone_of_node[n] == z))

    # DC power flow constraints
    fmap  = @constraint(model, [t in JH, l in JL], f[t, l] == sum(PTDF[l, n] * r[t, n] for n in JN))

    # thermal limits
    therm = @constraint(model, [t in JH, l in JL], -Fmax[l] <= f[t,l] <= Fmax[l])

    # system balance
    sbal  = @constraint(model, [t in JH], sum(r[t, n] for n in JN) == 0)

    # zonal capacity limits
    zonal_cap = @constraint(model, [t in JH, i in JI, z in JZ], g[t, i, z] <= AV[t, i, z] * (y[i, z] + sum(y_node[i, n] for n in JN if zone_of_node[n] == z)))

    # store metadata as a Dict — JuMP.Model expects model.ext to be a Dict-like object
    model.ext = Dict{Symbol,Any}(
        :sets => sets,
        :maps => maps,
        :params => params,
        :vars => Dict(
            :y => y, :y_bar => y_bar, :g => g, :g_bar => g_bar,
            :r => r, :f => f, :p => p, :d_inel => d_inel, :d_ela => d_ela, :ens => ens
        ),
        :constraint => Dict(
            :bal => bal, :agg => agg, :nbal => nbal, :cap => cap,
            :fmap => fmap, :therm => therm, :sbal => sbal,
            :zonal_cap => zonal_cap, :Ren_cap => Ren_cap,
            :alloc => alloc,
        ),
        :scalars => Dict(:ela => ela, :WTP => WTP),
        :expressions => Dict(:gen_cost => gen_cost, 
            :inv_cost => inv_cost, :neg_utility => neg_utility)
    )
    return model
end

function solve_and_save(
    data::Dict,
    load::DataFrame,
    pv::DataFrame,
    wind_on::DataFrame,
    nptdf::DataFrame,
    lines::DataFrame,
    weights; 
    output_dir = joinpath(@__DIR__, "Results", "Planner"),
)
    mkpath(output_dir)

    # Build and solve
    model = build_planner_eom(;
        data=data,
        load=load,
        pv=pv,
        wind_on=wind_on,
        nodal_ptdf_df=nptdf,
        lines=lines,
        weights=weights,
    )
    optimize!(model)

    term = termination_status(model)
    if term != MOI.OPTIMAL
        ps = primal_status(model)
        error("Solver did not return OPTIMAL. termination_status=$(term), primal_status=$(ps)")
    end

    ext    = model.ext
    sets   = ext[:sets]
    maps   = ext[:maps]
    params = ext[:params]
    expressions = ext[:expressions]

    JH, JI, JZ, JN, JL = sets.JH, sets.JI, sets.JZ, sets.JN, sets.JL
    zones, nodes, techs = sets.zones, sets.nodes, sets.techs

    vars = ext[:vars]
    @assert all(haskey(vars, s) for s in (:y, :y_bar, :g_bar, :p, :f, :d_inel, :d_ela)) "model.ext.vars is missing expected fields"

    y_var     = vars[:y]        # (I,Z)   installed zonal capacity
    ybar_var  = vars[:y_bar]    # (I,N)   nodal allocation of new capacity
    g_var     = vars[:g]        # (H,I,N) nodal generation
    gbar_var  = vars[:g_bar]    # (H,I,N) nodal dispatch
    p_var     = vars[:p]        # (H,Z)   zonal net positions
    f_var     = vars[:f]        # (H,L)   line flows
    dinel_var = vars[:d_inel]   # (H,Z)   inelastic demand served
    dela_var  = vars[:d_ela]    # (H,Z)   elastic demand served
    ens_var   = vars[:ens]      # (H,Z)   energy not served


    capacity     = value.(y_var)
    aux_capacity  = value.(ybar_var)
    gen          = value.(g_var)
    aux_dispatch = value.(gbar_var)
    net_pos = value.(p_var)
    flow     = value.(f_var)
    inel   = value.(dinel_var)
    elastic  = value.(dela_var)
    ens_var   = value.(ens_var)

    # Zonal prices (€/MWh) from balance duals, unweighted
    W = params.W
    rho = [dual(ext[:constraint][:bal][t, z]) / W[t,z] for t in JH, z in JZ]
    nbal_dual = [-dual(ext[:constraint][:nbal][t, n]) for t in JH, n in JN]


    # extract dual of alloc
    dual_alloc = [dual(ext[:constraint][:alloc][i, z]) for i in JI, z in JZ]
    println(DataFrame(dual_alloc, :auto))

    # Existing capacity C[i,z] from input data (for totals)
    Cz = zeros(length(JI), length(JZ))
    for (zname, zidx) in maps.zone_to_idx
        for (tname, pars) in data["Generators"][zname]
            i = maps.tech_to_idx[Symbol(tname)]
            Cz[i, zidx] = Float64(get_any(pars, "C"))
        end
    end



    for z in zones
        zidx = maps.zone_to_idx[z]
        df = DataFrame(Timestep = collect(JH))
        df.EOM_price = rho[:, zidx]

            # Get zonal generation per tech
            for tech in techs
                i = maps.tech_to_idx[tech]
                gen_t = gen[:, i, zidx]  # Direct access to zonal generation
                tstr = String(tech)
                df[!, Symbol("Gen_$(z)_$(tstr)")]                = gen_t
                df[!, Symbol("new_Capacity_Gen_$(z)_$(tstr)")]   = fill(capacity[i, zidx], length(JH))
                df[!, Symbol("Capacity_Gen_$(z)_$(tstr)")] = fill(Cz[i, zidx] + capacity[i, zidx], length(JH))
            end

        cons_total = inel[:, zidx] .+ elastic[:, zidx]
        # ens = value.(expressions[:ens])[:, zidx]

        df[!, Symbol("Cons_$(z)")]           = cons_total
        df[!, Symbol("Inelastic_Cons_$(z)")] = inel[:, zidx]
        df[!, Symbol("Elastic_Cons_$(z)")]   = elastic[:, zidx]
        df[!, :NetworkManager]               = -1 * net_pos[:, zidx]  # (imports positive)
        df[!, Symbol("ENS_Cons_$(z)")]       = ens_var[:, zidx]


        for col in names(df)
            if eltype(df[!, col]) <: Number
                df[!, col] = round.(df[!, col], digits=2)
            end
        end

        CSV.write(joinpath(output_dir, "planner_zone_$(z).csv"), df; delim=";")
    end

    price_df = DataFrame(Timestep = collect(JH))
    for (zidx, z) in enumerate(zones)
        price_df[!, Symbol("Price_$(z)")] = rho[:, zidx]
    end

    invest_df = DataFrame(Technology = String.(techs))
    for (zidx, z) in enumerate(zones)
        invest_df[!, Symbol("Invest_$(z)")] = capacity[:, zidx]
    end

    nodal_invest_df = DataFrame(
        Technology = repeat(String.(techs), sets.N),
        Node       = repeat(String.(String.(nodes)), inner = sets.I),
        Allocation = vec(aux_capacity),
    )

    dispatch_df = DataFrame(
        Timestep   = repeat(collect(JH), outer = sets.I * sets.N),
        Technology = repeat(repeat(String.(techs), inner = sets.H), outer = sets.N),
        Node       = repeat(String.(String.(nodes)), inner = sets.H * sets.I),
        Generation = vec(aux_dispatch),
    )


    flow_df = DataFrame(Timestep = repeat(collect(JH), outer = sets.L))
    line_ids = :line_id ∈ names(lines) ? lines[!, :line_id] : 1:nrow(lines)
    flow_df[!, :Line] = repeat(line_ids, inner = sets.H)
    flow_df[!, :Flow] = vec(flow)

    return Dict(
        :status               => term,
        :objective            => objective_value(model),
        :prices_df            => price_df,
        :investments_df       => invest_df,
        :nodal_investments_df => nodal_invest_df,
        :dispatch_df          => dispatch_df,
        :flows_df             => flow_df,
        :out_dir              => output_dir,
        :model                => model,
        :nbal_dual            => nbal_dual,  # Add these to the returned dictionary
        :rho                  => rho
    )
end

inputs = read_inputs(@__DIR__)
cp_eom = solve_and_save(inputs.data, inputs.load, inputs.pv, inputs.wind_on, inputs.nptdf, inputs.lines, inputs.weights; output_dir=joinpath(@__DIR__, "Results", "Planner"))
@show objective_value(cp_eom[:model])
# print(cp_eom)


JH = cp_eom[:model].ext[:sets].JH
JZ = cp_eom[:model].ext[:sets].JZ
JN = cp_eom[:model].ext[:sets].JN
zone_of_node = cp_eom[:model].ext[:maps].zone_of_node
nbal_dual = cp_eom[:nbal_dual]  
rho = cp_eom[:rho]          
y_nodal = cp_eom[:model].ext[:vars][:y_bar]

println(value.(y_nodal))


# for t in JH
#     for n in JN
#         z = zone_of_node[n]
#         nodal_price = nbal_dual[t, n]
#         zonal_price = rho[t, z]
#         # Print all nodal and zonal prices without condition
#         println("Time $t, Node $n (Zone $z): $(round(nodal_price, digits=2)), Zonal price = $(round(zonal_price, digits=2))")
#     end
# end
