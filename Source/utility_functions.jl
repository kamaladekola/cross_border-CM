# parse agent names
function parse_agent_name(agent_name::String)
    parts = split(agent_name, "_")
    if length(parts) == 3 && parts[1] == "Gen"
        return (parts[2], parts[3])  # (zone, tech)
    elseif length(parts) == 2 && parts[1] == "Cons"
        return (parts[2], nothing)  # (zone, nothing)
    elseif length(parts) == 3 && parts[1] == "IC"
        return (parts[2], parts[3]) # (z1, z2)
    else
        error("Agent name format not recognized: $agent_name")
    end
end

# check if agent is an interconnector
function is_interconnector(agent_name::String)::Bool
    parts = split(agent_name, "_")
    return length(parts) == 3 && parts[1] == "IC"
end

# get the sign of the interconnector flow 
function sign_for_zone(IC_m::String, z::String)
    z1, z2 = parse_agent_name(IC_m)  # "IC_A_B" => ("A","B")
    # -1 for export, +1 for import, 0 for no flow --> still some confusion with signs and convergence
    if z == z1
        return -1
    elseif z == z2
        return +1
    else
        return 0
    end
end
