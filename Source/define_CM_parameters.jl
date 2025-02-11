function define_CM_parameters!(CM::Dict, data::Dict, ts::DataFrame, scenario_overview_row::DataFrameRow)
    CM["local_price"] = data["CapacityMarket"]["local"]
    CM["foreign_price"] = data["CapacityMarket"]["foreign"]
    CM["demand_local"] = data["CapacityMarket"]["demand"]["A"]
    CM["demand_foreign"] = data["CapacityMarket"]["foreign_demand"]["A"]
    return CM
end
