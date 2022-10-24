# define glider data types for SeaExplorer glider
# gong@vims.edu: 2022-10-23
module seaexplorer_types

using Dates

mutable struct NAV_RT
    t::Array{DateTime};
    z::Array{AbstractFloat};
    lon::Array{AbstractFloat};
    lat::Array{AbstractFloat};
    NavState::Array{Int64};
    SecurityLevel::Array{Int64};
    Heading::Array{AbstractFloat};
    Declination::Array{Int64};
    Pitch::Array{AbstractFloat};
    Roll::Array{AbstractFloat};
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
    t::Array{DateTime};
    z::Array{AbstractFloat};
    lon::Array{AbstractFloat};
    lat::Array{AbstractFloat};
    nav_resource::Array{AbstractFloat};
    ad2cp_time::Array{DateTime};
    ad2cp_heading::Array{AbstractFloat};
    ad2cp_pitch::Array{AbstractFloat};
    ad2cp_roll::Array{AbstractFloat};
    ad2cp_pressure::Array{AbstractFloat};
    ad2cp_alt::Array{AbstractFloat};
    ad2cp_v1_cn1::Array{AbstractFloat};
    ad2cp_v2_cn1::Array{AbstractFloat};
    ad2cp_v3_cn1::Array{AbstractFloat};
    ad2cp_v4_cn1::Array{AbstractFloat};
    ad2cp_v1_cn2::Array{AbstractFloat};
    ad2cp_v2_cn2::Array{AbstractFloat};
    ad2cp_v3_cn2::Array{AbstractFloat};
    ad2cp_v4_cn2::Array{AbstractFloat};
    ad2cp_v1_cn3::Array{AbstractFloat};
    ad2cp_v2_cn3::Array{AbstractFloat};
    ad2cp_v3_cn3::Array{AbstractFloat};
    ad2cp_v4_cn3::Array{AbstractFloat};
    ad2cp_v1_cn4::Array{AbstractFloat};
    ad2cp_v2_cn4::Array{AbstractFloat};
    ad2cp_v3_cn4::Array{AbstractFloat};
    ad2cp_v4_cn4::Array{AbstractFloat};
    ad2cp_v1_cn5::Array{AbstractFloat};
    ad2cp_v2_cn5::Array{AbstractFloat};
    ad2cp_v3_cn5::Array{AbstractFloat};
    ad2cp_v4_cn5::Array{AbstractFloat};
    ad2cp_v1_cn6::Array{AbstractFloat};
    ad2cp_v2_cn6::Array{AbstractFloat};
    ad2cp_v3_cn6::Array{AbstractFloat};
    ad2cp_v4_cn6::Array{AbstractFloat};
    flbbcd_chl_count::Array{AbstractFloat};
    flbbcd_chl_scaled::Array{AbstractFloat};
    flbbcd_bb_700_count::Array{AbstractFloat};
    flbbcd_bb_700_scaled::Array{AbstractFloat};
    flbbcd_cdom_count::Array{AbstractFloat};
    flbbcd_cdom_scaled::Array{AbstractFloat};
    legato_conductivity::Array{AbstractFloat};
    legato_temperature::Array{AbstractFloat};
    legato_pressure::Array{AbstractFloat};
    legato_salinity::Array{AbstractFloat};
    legato_condtemp::Array{AbstractFloat};
    mr1000g_t1_avg::Array{AbstractFloat};
    mr1000g_t2_avg::Array{AbstractFloat};
    mr1000g_sh1_std::Array{AbstractFloat};
    mr1000g_sh2_std::Array{AbstractFloat};
    mr1000g_press_avg::Array{AbstractFloat};
    mr1000g_incly_avg::Array{AbstractFloat};
    mr1000g_eps1::Array{AbstractFloat};
    mr1000g_qc1::Array{AbstractFloat};
    mr1000g_eps2::Array{AbstractFloat};
    mr1000g_qc2::Array{AbstractFloat};
end

mutable struct NAV
    t::Array{DateTime};
    z::Array{AbstractFloat};
    lon::Array{AbstractFloat};
    lat::Array{AbstractFloat};
end

mutable struct LEGATO
    t::Array{DateTime};
    p::Array{AbstractFloat};
    temp::Array{AbstractFloat};
    cond::Array{AbstractFloat};
    condtemp::Array{AbstractFloat};
    salt::Array{AbstractFloat};
end

mutable struct FLBBCD
    t::Array{DateTime};
    chla::Array{AbstractFloat};
    cdom::Array{AbstractFloat};
    bb700::Array{AbstractFloat};
end

mutable struct AD2CP
    t::Array{DateTime};
    heading::Array{AbstractFloat};
    pitch::Array{AbstractFloat};
    roll::Array{AbstractFloat};
    p::Array{AbstractFloat};
    alt::Array{AbstractFloat};
    v1_ch1::Array{AbstractFloat};
    v2_ch1::Array{AbstractFloat};
    v3_ch1::Array{AbstractFloat};
    v4_ch1::Array{AbstractFloat};
    v1_ch2::Array{AbstractFloat};
    v2_ch2::Array{AbstractFloat};
    v3_ch2::Array{AbstractFloat};
    v4_ch2::Array{AbstractFloat};
    v1_ch3::Array{AbstractFloat};
    v2_ch3::Array{AbstractFloat};
    v3_ch3::Array{AbstractFloat};
    v4_ch3::Array{AbstractFloat};
    v1_ch4::Array{AbstractFloat};
    v2_ch4::Array{AbstractFloat};
    v3_ch4::Array{AbstractFloat};
    v4_ch4::Array{AbstractFloat};
    v1_ch5::Array{AbstractFloat};
    v2_ch5::Array{AbstractFloat};
    v3_ch5::Array{AbstractFloat};
    v4_ch5::Array{AbstractFloat};
    v1_ch6::Array{AbstractFloat};
    v2_ch6::Array{AbstractFloat};
    v3_ch6::Array{AbstractFloat};
    v4_ch6::Array{AbstractFloat};
end

mutable struct MR1000G
    t::Array{DateTime};
    mr1000g_t1_avg::Array{AbstractFloat};
    mr1000g_t2_avg::Array{AbstractFloat};
    mr1000g_sh1_std::Array{AbstractFloat};
    mr1000g_sh2_std::Array{AbstractFloat};
    mr1000g_press_avg::Array{AbstractFloat};
    mr1000g_incly_avg::Array{AbstractFloat};
    mr1000g_eps1::Array{AbstractFloat};
    mr1000g_qc1::Array{AbstractFloat};
end

mutable struct SeaExplorer
    nav::Array{NAV};
    ctd::Array{LEGATO};
    flbbcd::Array{FLBBCD};
    ad2cp::Array{AD2CP}
    nav1d::NAV;
    ctd1d::LEGATO;
    flbbcd1d::FLBBCD
    ad2cp1d::AD2CP;
end

mutable struct SeaExplorerRT
    nav::Array{NAV_RT};
    pld::Array{PLD_RT};
    nav1d::Array{NAV_RT};
    pld1d::Array{PLD_RT};
end

end