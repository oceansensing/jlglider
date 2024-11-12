module slocumType

mutable struct engStruct
    glidertype::String
    gliderSN::Int
    glidername::String
    missionID::Int
    project::String
    t::Array{Float64}
    p::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
end

mutable struct ctdStruct
    glidertype::String
    gliderSN::Int
    glidername::String
    missionID::Int
    project::String
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
    glidertype::String
    gliderSN::Int
    glidername::String
    missionID::Int
    project::String
    t::Array{Float64}
    p::Array{Float64}
    z::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
    var::Array{Float64}
end  

end