# define glider data types for SeaExplorer glider
# gong@vims.edu: 2024-10-28 this is a work in progress trying to merge Slocum and SeaExplorer glider types, not usable yet, type logic doesn't even make sense, just a start.
module gliderType

using Dates

abstract type PLATFORM end
abstract type INSTRUMENT end
abstract type GLIDER <: PLATFORM end
abstract type CTD <: INSTRUMENT end
abstract type SEAEXPLORER <: GLIDER & INSTRUMENT & NAVIGATION end

mutable struct GliderMission
    project::String
    gliderType::String
    gliderID::Int64
    gliderName::String
    startDate::Date
end

mutable struct NAVIGATION 
    t::Array{DateTime};
    z::Array{Float64};
    lon::Array{Float64};
    lat::Array{Float64};
end

struct SEAEXPLORER <: GLIDER
    mission::GliderMission
    nav::NAVIGATION
end

struct SLOCUM <: GLIDER

end

mutable struct LEGATO end
mutable struct SEABIRD end


mutable struct NAV_RT
    gliderSN::Int64
    missionID::Int64
    project::String    
    yo::Array{Int64};
    t::Array{DateTime};
    z::Array{Float64};
    lon::Array{Float64};
    lat::Array{Float64};
    NavState::Array{Int64};
    SecurityLevel::Array{Int64};
    Heading::Array{Float64};
    Declination::Array{Int64};
    Pitch::Array{Float64};
    Roll::Array{Float64};
    DeadReckoning::Array{Int64};
    DesiredH::Array{Int64};
    BallastCmd::Array{Int64};
    BallastPos::Array{Float64};
    LinCmd::Array{Float64};
    LinPos::Array{Float64};
    AngCmd::Array{Int64};
    AngPos::Array{Float64};
    Voltage::Array{Float64};
    Altitude::Array{Float64};
end

mutable struct PLD_RT
    gliderSN::Int64
    missionID::Int64
    project::String    
    yo::Array{Int64};
    t::Array{DateTime};
    z::Array{Float64};
    lon::Array{Float64};
    lat::Array{Float64};
    nav_resource::Array{Float64};
    ad2cp_time::Array{DateTime};
    ad2cp_heading::Array{Float64};
    ad2cp_pitch::Array{Float64};
    ad2cp_roll::Array{Float64};
    ad2cp_pressure::Array{Float64};
    ad2cp_alt::Array{Float64};
    ad2cp_v1_cn1::Array{Float64};
    ad2cp_v2_cn1::Array{Float64};
    ad2cp_v3_cn1::Array{Float64};
    ad2cp_v4_cn1::Array{Float64};
    ad2cp_v1_cn2::Array{Float64};
    ad2cp_v2_cn2::Array{Float64};
    ad2cp_v3_cn2::Array{Float64};
    ad2cp_v4_cn2::Array{Float64};
    ad2cp_v1_cn3::Array{Float64};
    ad2cp_v2_cn3::Array{Float64};
    ad2cp_v3_cn3::Array{Float64};
    ad2cp_v4_cn3::Array{Float64};
    ad2cp_v1_cn4::Array{Float64};
    ad2cp_v2_cn4::Array{Float64};
    ad2cp_v3_cn4::Array{Float64};
    ad2cp_v4_cn4::Array{Float64};
    ad2cp_v1_cn5::Array{Float64};
    ad2cp_v2_cn5::Array{Float64};
    ad2cp_v3_cn5::Array{Float64};
    ad2cp_v4_cn5::Array{Float64};
    ad2cp_v1_cn6::Array{Float64};
    ad2cp_v2_cn6::Array{Float64};
    ad2cp_v3_cn6::Array{Float64};
    ad2cp_v4_cn6::Array{Float64};
    flbbcd_chl_count::Array{Float64};
    flbbcd_chl_scaled::Array{Float64};
    flbbcd_bb_700_count::Array{Float64};
    flbbcd_bb_700_scaled::Array{Float64};
    flbbcd_cdom_count::Array{Float64};
    flbbcd_cdom_scaled::Array{Float64};
    legato_conductivity::Array{Float64};
    legato_temperature::Array{Float64};
    legato_pressure::Array{Float64};
    legato_salinity::Array{Float64};
    legato_condtemp::Array{Float64};
    mr1000g_t1_avg::Array{Float64};
    mr1000g_t2_avg::Array{Float64};
    mr1000g_sh1_std::Array{Float64};
    mr1000g_sh2_std::Array{Float64};
    mr1000g_press_avg::Array{Float64};
    mr1000g_incly_avg::Array{Float64};
    mr1000g_eps1::Array{Float64};
    mr1000g_qc1::Array{Float64};
    mr1000g_eps2::Array{Float64};
    mr1000g_qc2::Array{Float64};
    ad2cp_Unorth::Array{Float64};
    ad2cp_Ueast::Array{Float64};
    ad2cp_Utot::Array{Float64};
    ad2cp_Udir::Array{Float64};
    ad2cp_qf::Array{Float64};
end

mutable struct NAV
    t::Array{DateTime};
    z::Array{Float64};
    lon::Array{Float64};
    lat::Array{Float64};
end

mutable struct LEGATO
    t::Array{DateTime};
    p::Array{Float64};
    temp::Array{Float64};
    cond::Array{Float64};
    condtemp::Array{Float64};
    salt::Array{Float64};
end

mutable struct FLBBCD
    t::Array{DateTime};
    chla::Array{Float64};
    cdom::Array{Float64};
    bb700::Array{Float64};
end

mutable struct AD2CP
    t::Array{DateTime};
    heading::Array{Float64};
    pitch::Array{Float64};
    roll::Array{Float64};
    p::Array{Float64};
    alt::Array{Float64};
    v1_ch1::Array{Float64};
    v2_ch1::Array{Float64};
    v3_ch1::Array{Float64};
    v4_ch1::Array{Float64};
    v1_ch2::Array{Float64};
    v2_ch2::Array{Float64};
    v3_ch2::Array{Float64};
    v4_ch2::Array{Float64};
    v1_ch3::Array{Float64};
    v2_ch3::Array{Float64};
    v3_ch3::Array{Float64};
    v4_ch3::Array{Float64};
    v1_ch4::Array{Float64};
    v2_ch4::Array{Float64};
    v3_ch4::Array{Float64};
    v4_ch4::Array{Float64};
    v1_ch5::Array{Float64};
    v2_ch5::Array{Float64};
    v3_ch5::Array{Float64};
    v4_ch5::Array{Float64};
    v1_ch6::Array{Float64};
    v2_ch6::Array{Float64};
    v3_ch6::Array{Float64};
    v4_ch6::Array{Float64};
end

mutable struct MR1000G
    t::Array{DateTime};
    mr1000g_t1_avg::Array{Float64};
    mr1000g_t2_avg::Array{Float64};
    mr1000g_sh1_std::Array{Float64};
    mr1000g_sh2_std::Array{Float64};
    mr1000g_press_avg::Array{Float64};
    mr1000g_incly_avg::Array{Float64};
    mr1000g_eps1::Array{Float64};
    mr1000g_qc1::Array{Float64};
end

mutable struct SeaExplorerTest
    nav::Array{NAV};
    ctd::Array{LEGATO};
    flbbcd::Array{FLBBCD};
    ad2cp::Array{AD2CP}
    nav1d::NAV;
    ctd1d::LEGATO;
    flbbcd1d::FLBBCD
    ad2cp1d::AD2CP;
end

mutable struct SeaExplorer
    nav::Array{NAV_RT};
    pld::Array{PLD_RT};
    nav1d::Array{NAV_RT};
    pld1d::Array{PLD_RT};
end

mutable struct SeaExplorerData
    gliderSN::Int64
    missionID::Int64
    project::String
    yo::Array{Int64}
    ns::Array{Int64}
    t::Array{DateTime}
    lon::Array{Float64}
    lat::Array{Float64}
    gpst::Array{DateTime}
    gpslon::Array{Float64}
    gpslat::Array{Float64}
    p::Array{Float64}
    z::Array{Float64}
    temp::Array{Float64}
    salt::Array{Float64}
    saltA::Array{Float64}
    ctemp::Array{Float64}
    sigma0::Array{Float64}
    spice0::Array{Float64}
    sndspd::Array{Float64}
 
    mr_eps1::Array{Float64}
    mr_eps2::Array{Float64}
    mr_qc1::Array{Float64}
    mr_qc2::Array{Float64}
    mr_sh1_std::Array{Float64}
    mr_sh2_std::Array{Float64}
    mr_t1_avg::Array{Float64}
    mr_t2_avg::Array{Float64}
    ad2cp_Ueast::Array{Float64}
    ad2cp_Unorth::Array{Float64}
    ad2cp_Utot::Array{Float64}
    ad2cp_Udir::Array{Float64}
    ad2cp_qf::Array{Float64}
    chla::Array{Float64}
    bb700::Array{Float64}
    cdom::Array{Float64}
    n2::Array{Float64}
    pmid::Array{Float64}
    zmid::Array{Float64}
    tmid::Array{DateTime}
end

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
    mission::String
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

end