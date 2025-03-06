# parse agent names
function parse_agent_name(agent_name::String)
    parts = split(agent_name, "_")
    if agent_name == "NetworkManager"
        return ("NetworkManager", nothing)
    elseif length(parts) == 3 && parts[1] == "Gen"
        return (parts[2], parts[3])  # (zone, tech)
    elseif length(parts) == 2 && parts[1] == "Cons"
        return (parts[2], nothing)  # (zone, nothing)
    else
        error("Agent name format not recognized: $agent_name")
    end
end


function get_agent_name(m::String)
    agent_name = split(m, "_")
    return agent_name[end]
end

