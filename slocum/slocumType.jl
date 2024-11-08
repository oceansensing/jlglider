module slocumType

mutable struct plotSetting
    #figoutdir::String
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
    mission::String
    glidername::String
    #resolution::Tuple{Int64, Int64}
    #markersize::Float64
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

mutable struct engStruct
    t::Array{Float64}
    p::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
end

mutable struct ctdStruct
    project::String
    gliderSN::Int
    glidername::String
    t::Array{Float64}
    p::Array{Float64}
    z::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
    temp::Array{Float64}
    cond::Array{Float64}
    salt::Array{Float64}
    ctemp::Array{Float64}
    saltA::Array{Float64}
    sigma0::Array{Float64}
    spice0::Array{Float64}
    sndspd::Array{Float64}
    castnum::Int
    qflag::Int
end

mutable struct sciStruct
    mission::String
    glidername::String    
    t::Array{Float64}
    p::Array{Float64}
    z::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
    var::Array{Float64}
end  

end