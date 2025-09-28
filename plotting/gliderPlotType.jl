module gliderPlotType

export plotSetting, plotStruct

mutable struct plotSetting
    pint::Float64
    iday::Float64
    ms::Float64
    tsms::Float64
    pres::Tuple{Int64, Int64}
    tspres::Tuple{Int64, Int64}
    fs::Float64
end

mutable struct plotStruct
    figoutdir::String
    project::String
    glidername::String
    lonmin::Float64
    lonmax::Float64
    latmin::Float64
    latmax::Float64
    pmin::Float64
    pmax::Float64
    zmin::Float64
    zmax::Float64
    tempmin::Float64
    tempmax::Float64
    condmin::Float64
    condmax::Float64
    saltmin::Float64
    saltmax::Float64
    sigma0min::Float64
    sigma0max::Float64
    spice0min::Float64
    spice0max::Float64
    sndspdmin::Float64
    sndspdmax::Float64
end

end