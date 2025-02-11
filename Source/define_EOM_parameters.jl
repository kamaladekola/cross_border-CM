function define_EOM_parameters!(EOM::Dict, data::Dict, load::DataFrame, scenario_overview_row::DataFrameRow, zones::Vector{String})

    # EOM["D"] = Dict(z => load[!, Symbol("LOAD_", z)][1:data["General"]["nTimesteps"]] for z in zones) # Load timeseries (MWh)

    return EOM
end

# not used