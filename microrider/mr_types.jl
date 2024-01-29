module MR_types

mutable struct MicroRiderRaw
    fullPath::String
    fs_fast::AbstractFloat
    fs_slow::AbstractFloat
    header_version::Float64
    t_slow::Array{Float64}
    t_fast::Array{Float64}
    setupfilestr::String
    cfgobj::Dict{String, Any}
    header::Array{Float64}
    filetime::Float64
    date::String
    time::String
    Gnd::Array{Float64}
    Ax::Array{Float64}
    Ay::Array{Float64}
    T1::Array{Float64}
    T1_dT1::Array{Float64}
    T2::Array{Float64}
    T2_dT2::Array{Float64}
    sh1::Array{Float64}
    sh2::Array{Float64}
    P::Array{Float64}
    P_dP::Array{Float64}
    PV::Array{Float64}
    V_Bat::Array{Float64}
    Incl_Y::Array{Float64}
    Incl_X::Array{Float64}
    Incl_T::Array{Float64}
    odas_version::Float64
    vehicle_info::Dict{String, Any}
    t_fast_YD::Array{Float64}
    t_slow_YD::Array{Float64}
    Year::Float64
    Month::Float64
    Day::Float64
    Hour::Float64
    Minute::Float64
    Second::Float64
    Milli::Float64
    T1_slow::Array{Float64}
    T1_fast::Array{Float64}
    T2_slow::Array{Float64}
    T2_fast::Array{Float64}
    P_slow::Array{Float64}
    P_fast::Array{Float64}
    temperature_fast::Array{Float64}
    W_slow::Array{Float64}
    W_fast::Array{Float64}
    speed_slow::Array{Float64}
    speed_fast::Array{Float64}
    gradT1::Array{Float64}
    gradT2::Array{Float64}
    input_parameters::Dict{String, Any}
    params::Dict{String, Any}
end

mutable struct MicroRider{T}
    mr::T
    z::Array{Float64}
end

end