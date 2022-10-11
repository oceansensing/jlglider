# This script load SeaExplorer data file
# 2022-07-16: Donglai Gong
# 2022-10-10: DG adding loading of NAV and RT data, reading gzipped files

using Glob, DataFrames, CSV, Dates, Missings
import TranscodingStreams, CodecZlib
import GZip

# creating glider data types
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
    ad2cp_time::Array{AbstractFloat};
    ad2cp_heading::Array{AbstractFloat};
    ad2cp_pitch::Array{AbstractFloat};
    ad2cp_roll::Array{AbstractFloat};
    ad2cp_pressure::Array{AbstractFloat};
    ad2cp_alt::Array{AbstractFloat};
    ad2cp_v1_ch1::Array{AbstractFloat};
    ad2cp_v2_ch1::Array{AbstractFloat};
    ad2cp_v3_ch1::Array{AbstractFloat};
    ad2cp_v4_ch1::Array{AbstractFloat};
    ad2cp_v1_ch2::Array{AbstractFloat};
    ad2cp_v2_ch2::Array{AbstractFloat};
    ad2cp_v3_ch2::Array{AbstractFloat};
    ad2cp_v4_ch2::Array{AbstractFloat};
    ad2cp_v1_ch3::Array{AbstractFloat};
    ad2cp_v2_ch3::Array{AbstractFloat};
    ad2cp_v3_ch3::Array{AbstractFloat};
    ad2cp_v4_ch3::Array{AbstractFloat};
    ad2cp_v1_ch4::Array{AbstractFloat};
    ad2cp_v2_ch4::Array{AbstractFloat};
    ad2cp_v3_ch4::Array{AbstractFloat};
    ad2cp_v4_ch4::Array{AbstractFloat};
    ad2cp_v1_ch5::Array{AbstractFloat};
    ad2cp_v2_ch5::Array{AbstractFloat};
    ad2cp_v3_ch5::Array{AbstractFloat};
    ad2cp_v4_ch5::Array{AbstractFloat};
    ad2cp_v1_ch6::Array{AbstractFloat};
    ad2cp_v2_ch6::Array{AbstractFloat};
    ad2cp_v3_ch6::Array{AbstractFloat};
    ad2cp_v4_ch6::Array{AbstractFloat};
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

# define a function that converts missings in an array to NaN
function missing2nan(varin)
    varout = Float64.(collect(Missings.replace(varin, NaN)));
end

function load_SEAdata_rt(glidername::String, mission::String, navdir::String, scidir::String)
    nav_rt = NAV_RT[];
    pld_rt = PLD_RT[];

    missionroot = glidername * "." * mission;
    gliroot_rt = missionroot * "." * "gli.sub.";
    pldroot_rt = missionroot * "." * "pld1.sub.";

    glilist_rt = Glob.glob(gliroot_rt * "*", navdir);
    pldlist_rt = Glob.glob(pldroot_rt * "*", scidir);

    # time data format
    timeformat = "dd/mm/yyyy HH:MM:SS.sss"
    
    # initiate NAV_RT 1-D variables
    t1d = [];
    z1d = [];
    lon1d = [];
    lat1d = [];
    NavState1d = [];
    SecurityLevel1d = [];
    Heading1d = [];
    Declination1d = [];
    Pitch1d = [];
    Roll1d = [];
    DeadReckoning1d = [];
    DesiredH1d = [];
    BallastCmd1d = [];
    BallastPos1d = [];
    LinCmd1d = [];
    LinPos1d = [];
    AngCmd1d = [];
    AngPos1d = [];
    Voltage1d = [];
    Altitude1d = [];    

    # loop through the list of realtime navigation (nav) data files
    for i = 1:length(glilist_rt)
        display(i)
        navfilepath = navdir * gliroot_rt * string(i) * ".gz"; 
        print(navfilepath * "\n")
        df = CSV.read(navfilepath, header=1, delim=";", DataFrame, buffer_in_memory=true);

        t = DateTime.(df.Timestamp, timeformat);
        z = df.Depth;
        navlon = df.Lon;
        navlat = df.Lat;
        lon = trunc.(navlon ./ 100) + (navlon .- trunc.(navlon ./ 100)*100) / 60;
        lat = trunc.(navlat ./ 100) + (navlat .- trunc.(navlat ./ 100)*100) / 60;

        NavState = df.NavState;
        SecurityLevel = df.SecurityLevel;
        Heading = df.Heading;
        Declination = df.Declination;
        Pitch = df.Pitch;
        Roll = df.Roll;
        DeadReckoning = df.DeadReckoning;
        DesiredH = df.DesiredH;
        BallastCmd = df.BallastCmd;
        BallastPos = df.BallastPos;
        LinCmd = df.LinCmd;
        LinPos = df.LinPos;
        AngCmd = df.AngCmd;
        AngPos = df.AngPos;
        Voltage = df.Voltage;
        Altitude = df.Altitude;

        gt1d = cat(gt1d, t, dims=1);
        gz1d = cat(gz1d, z, dims=1);
        glon1d = cat(glon1d, lon, dims=1);
        glat1d = cat(glat1d, lat, dims=1);
        NavState1d = cat(NavState1d, NavState, dims=1);
        SecurityLevel1d = cat(SecurityLevel1d, SecurityLevel, dims=1);
        Heading1d = cat(Heading1d, Heading, dims=1);
        Declination1d = cat(Declination1d, Declination, dims=1);
        Pitch1d = cat(Pitch1d, Pitch, dims=1);
        Roll1d = cat(Roll1d, Roll, dims=1);
        DeadReckoning1d = cat(DeadReckoning1d, DeadReckoning, dims=1);
        DesiredH1d = cat(DesiredH1d, DesiredH, dims=1);
        BallastCmd1d = cat(BallastCmd1d, BallastCmd, dims=1);
        BallastPos1d = cat(BallastPos1d, BallastPos, dims=1);
        LinCmd1d = cat(LinCmd1d, LinCmd, dims=1);
        LinPos1d = cat(LinPos1d, LinPos, dims=1);
        AngCmd1d = cat(AngCmd1d, AngCmd, dims=1);
        AngPos1d = cat(AngPos1d, AngPos, dims=1);
        Voltage1d = cat(Voltage1d, Voltage, dims=1);
        Altitude1d = cat(Altitude1d, Altitude, dims=1);

        push!(nav_rt, NAV_RT(t, z, lon, lat, NavState, SecurityLevel, Heading, Declination, Pitch, Roll, DeadReckoning, DesiredH, BallastCmd, BallastPos, LinCmd, LinPos, AngCmd, AngPos, Voltage, Altitude));
    end
    nav1d_rt = NAV_RT(t1d, z1d, lon1d, lat1d, NavState1d, SecurityLevel1d, Heading1d, Declination1d, Pitch1d, Roll1d, DeadReckoning1d, DesiredH1d, BallastCmd1d, BallastPos1d, LinCmd1d, LinPos1d, AngCmd1d, AngPos1d, Voltage1d, Altitude1d);


    # initiate PLD_RT 1-D variables
    t1d = [];
    z1d = [];
    lon1d = [];
    lat1d = [];
    nav_resource1d = [];
    ad2cp_time1d = [];
    ad2cp_heading1d = [];
    ad2cp_pitch1d = [];
    ad2cp_roll1d = [];
    ad2cp_pressure1d = [];
    ad2cp_alt1d = [];
    ad2cp_v1_ch1_1d = [];
    ad2cp_v2_ch1_1d = [];
    ad2cp_v3_ch1_1d = [];
    ad2cp_v4_ch1_1d = [];
    ad2cp_v1_ch2_1d = [];
    ad2cp_v2_ch2_1d = [];
    ad2cp_v3_ch2_1d = [];
    ad2cp_v4_ch2_1d = [];
    ad2cp_v1_ch3_1d = [];
    ad2cp_v2_ch3_1d = [];
    ad2cp_v3_ch3_1d = [];
    ad2cp_v4_ch3_1d = [];
    ad2cp_v1_ch4_1d = [];
    ad2cp_v2_ch4_1d = [];
    ad2cp_v3_ch4_1d = [];
    ad2cp_v4_ch4_1d = [];
    ad2cp_v1_ch5_1d = [];
    ad2cp_v2_ch5_1d = [];
    ad2cp_v3_ch5_1d = [];
    ad2cp_v4_ch5_1d = [];
    ad2cp_v1_ch6_1d = [];
    ad2cp_v2_ch6_1d = [];
    ad2cp_v3_ch6_1d = [];
    ad2cp_v4_ch6_1d = [];
    flbbcd_chl_count_1d = [];
    flbbcd_chl_scaled_1d = [];
    flbbcd_bb_700_count_1d =[];
    flbbcd_bb_700_scaled_1d = [];
    flbbcd_cdom_count_1d = [];
    flbbcd_cdom_scaled_1d = [];
    legato_conductivity_1d = [];
    legato_temperature_1d = [];
    legato_pressure_1d = [];
    legato_salinity_1d = [];
    legato_condtemp_1d = [];
    mr1000g_t1_avg_1d = [];
    mr1000g_t2_avg_1d = [];
    mr1000g_sh1_std_1d = [];
    mr1000g_sh2_std_1d = [];
    mr1000g_press_avg_1d = [];
    mr1000g_incly_avg_1d = [];
    mr1000g_eps1_1d = [];
    mr1000g_qc1_1d = [];
    mr1000g_eps2_1d = [];
    mr1000g_qc2_1d = [];

    # loop through the list of realtime payload (pld) data files
    for i = 1:length(pldlist_rt)
        display(i)
        pldfilepath = scidir * pldroot_rt * string(i) * ".gz"; 
        print(pldfilepath * "\n")
        df = CSV.read(pldfilepath, header=1, delim=";", DataFrame, buffer_in_memory=true);

        # extract location data from data frame
        t = DateTime.(df.PLD_REALTIMECLOCK, timeformat);
        navlon = df.NAV_LONGITUDE;
        navlat = df.NAV_LATITUDE;
        lon = trunc.(navlon ./ 100) + (navlon .- trunc.(navlon ./ 100)*100) / 60;
        lat = trunc.(navlat ./ 100) + (navlat .- trunc.(navlat ./ 100)*100) / 60;
        lon = missing2nan(lon);
        lat = missing2nan(lat);
        z = missing2nan(df.NAV_DEPTH);
        
        push!(pld_rt, PLD_RT(t, z, lon, lat, nav_resource, ad2cp_time, ad2cp_heading, ad2cp_pitch, ad2cp_roll, ad2cp_pressure, ad2cp_alt, ad2cp_v1_ch1, ad2cp_v2_ch1, ad2cp_v3_ch1, ad2cp_v4_ch1, ad2cp_v1_ch2, ad2cp_v2_ch2, ad2cp_v3_ch2, ad2cp_v4_ch2, ad2cp_v1_ch3, ad2cp_v2_ch3, ad2cp_v3_ch3, ad2cp_v4_ch3, ad2cp_v1_ch4, ad2cp_v2_ch4, ad2cp_v3_ch4, ad2cp_v4_ch4, ad2cp_v1_ch5, ad2cp_v2_ch5, ad2cp_v3_ch5, ad2cp_v4_ch5, ad2cp_v1_ch6, ad2cp_v2_ch6, ad2cp_v3_ch6, ad2cp_v4_ch6, flbbcd_chl_count, flbbcd_chl_scaled, flbbcd_bb_700_count, flbbcd_bb_700_scaled, flbbcd_cdom_count, flbbcd_cdom_scaled, legato_conductivity, legato_temperature, legato_pressure, legato_salinity, legato_condtemp, mr1000g_t1_avg, mr1000g_t2_avg, mr1000g_sh1_std, mr1000g_sh2_std, mr1000g_press_avg, mr1000g_incly_avg, mr1000g_eps1, mr1000g_qc1, mr1000g_eps2, mr1000g_qc2));    
    end
    pld1d_rt = PLD_RT(t1d, z1d, lon1d, lat1d, nav_resource1d, ad2cp_time1d, ad2cp_heading1d, ad2cp_pitch1d, ad2cp_roll1d, ad2cp_pressure1d, ad2cp_alt1d, ad2cp_v1_ch1_1d, ad2cp_v2_ch1_1d, ad2cp_v3_ch1_1d, ad2cp_v4_ch1_1d, ad2cp_v1_ch2_1d, ad2cp_v2_ch2_1d, ad2cp_v3_ch2_1d, ad2cp_v4_ch2_1d, ad2cp_v1_ch3_1d, ad2cp_v2_ch3_1d, ad2cp_v3_ch3_1d, ad2cp_v4_ch3_1d, ad2cp_v1_ch4_1d, ad2cp_v2_ch4_1d, ad2cp_v3_ch4_1d, ad2cp_v4_ch4_1d, ad2cp_v1_ch5_1d, ad2cp_v2_ch5_1d, ad2cp_v3_ch5_1d, ad2cp_v4_ch5_1d, ad2cp_v1_ch6_1d, ad2cp_v2_ch6_1d, ad2cp_v3_ch6_1d, ad2cp_v4_ch6_1d, flbbcd_chl_count_1d, flbbcd_chl_scaled_1d, flbbcd_bb_700_count_1, flbbcd_bb_700_scaled_1d, flbbcd_cdom_count_1d, flbbcd_cdom_scaled_1d, legato_conductivity_1d, legato_temperature_1d, legato_pressure_1d, legato_salinity_1d, legato_condtemp_1d, mr1000g_t1_avg_1d, mr1000g_t2_avg_1d, mr1000g_sh1_std_1d, mr1000g_sh2_std_1d, mr1000g_press_avg_1d, mr1000g_incly_avg_1d, mr1000g_eps1_1d, mr1000g_qc1_1d, mr1000g_eps2_1d, mr1000g_qc2_1d);

    # combinating NAV_RT and PLD_RT data into one glider data structure
    gliderRT = SeaExplorerRT(nav_rt, pld_rt, nav1d_rt, pld1d_rt);

    return gliderRT
end

function load_SEAdata_raw(glidername::String, mission::String, navdir::String, scidir::String)

    # initilized data strucctures
    nav = NAV[];
    ctd = LEGATO[];
    flbbcd = FLBBCD[];
    ad2cp = AD2CP[];
    #glider = SeaExplorer[];

    # time data format
    timeformat = "dd/mm/yyyy HH:MM:SS.sss"

    # set directory paths
    missionroot = glidername * "." * mission;
    pldroot_raw = missionroot * "." * "pld1.raw.";
    ad2cproot_raw = missionroot * "." * "ad2cp.raw.";
    legatoroot_raw = missionroot * "." * "legato.raw.";

    # load data file lists
    pldlist_raw = Glob.glob(pldroot_raw * "*", scidir);
    ad2cplist_raw = Glob.glob(ad2cproot_raw * "*", scidir);
    legatolist_raw = Glob.glob(legatoroot_raw * "*", scidir);

    # initiate NAV 1-D variables 
    t1d = [];
    lon1d = [];
    lat1d = [];
    z1d = [];

    # initiate Legato 1-D variables
    temp1d = [];
    condtemp1d = [];
    cond1d = [];
    salt1d = [];
    p1d = [];

    # initiate FLBBCD 1-D variables
    chla1d = [];
    cdom1d = [];
    bb1d = [];

    # initiate ADCP 1-d variables
    ad2cp_time1d = [];
    ad2cp_heading1d = [];
    ad2cp_pitch1d = [];
    ad2cp_roll1d = [];
    ad2cp_pressure1d = [];
    ad2cp_alt1d = [];
    ad2cp_v1_ch1_1d = [];
    ad2cp_v2_ch1_1d = [];
    ad2cp_v3_ch1_1d = [];
    ad2cp_v4_ch1_1d = [];
    ad2cp_v1_ch2_1d = [];
    ad2cp_v2_ch2_1d = [];
    ad2cp_v3_ch2_1d = [];
    ad2cp_v4_ch2_1d = [];
    ad2cp_v1_ch3_1d = [];
    ad2cp_v2_ch3_1d = [];
    ad2cp_v3_ch3_1d = [];
    ad2cp_v4_ch3_1d = [];
    ad2cp_v1_ch4_1d = [];
    ad2cp_v2_ch4_1d = [];
    ad2cp_v3_ch4_1d = [];
    ad2cp_v4_ch4_1d = [];
    ad2cp_v1_ch5_1d = [];
    ad2cp_v2_ch5_1d = [];
    ad2cp_v3_ch5_1d = [];
    ad2cp_v4_ch5_1d = [];
    ad2cp_v1_ch6_1d = [];
    ad2cp_v2_ch6_1d = [];
    ad2cp_v3_ch6_1d = [];
    ad2cp_v4_ch6_1d = [];    

    # loop through the list of payload (pld) data files
    for i = 1:length(pldlist_raw)
        display(i)
        print(scidir * pldroot_raw * string(i) * "\n")
        df = CSV.read(scidir * pldroot_raw * string(i), header=1, delim=";", DataFrame);

        # extract location data from data frame
        t = DateTime.(df.PLD_REALTIMECLOCK, timeformat);
        navlon = df.NAV_LONGITUDE;
        navlat = df.NAV_LATITUDE;
        lon = trunc.(navlon ./ 100) + (navlon .- trunc.(navlon ./ 100)*100) / 60;
        lat = trunc.(navlat ./ 100) + (navlat .- trunc.(navlat ./ 100)*100) / 60;
        lon = missing2nan(lon);
        lat = missing2nan(lat);
        z = missing2nan(df.NAV_DEPTH);

        # concatinate NAV data into 1D arrays
        t1d = cat(t1d, t, dims=1);
        lon1d = cat(lon1d, lon, dims=1);
        lat1d = cat(lat1d, lat, dims=1);
        z1d = cat(z1d, z, dims=1);

        # change missings in LEGATO data to NaN
        p = missing2nan(df.LEGATO_PRESSURE);
        temp = missing2nan(df.LEGATO_TEMPERATURE);
        condtemp = missing2nan(df.LEGATO_CONDTEMP);
        cond = missing2nan(df.LEGATO_CONDUCTIVITY);
        salt = missing2nan(df.LEGATO_SALINITY);

        # concatinate LEGATO data into 1D arrays
        p1d = cat(p1d, p, dims=1);
        temp1d = cat(temp1d, temp, dims=1);
        condtemp1d = cat(condtemp1d, condtemp, dims=1);
        cond1d = cat(cond1d, cond, dims=1);
        salt1d = cat(salt1d, salt, dims=1);

        # change missings in FLBBCD to NaN
        chla = missing2nan(df.FLBBCD_CHL_SCALED);
        cdom = missing2nan(df.FLBBCD_CDOM_SCALED);
        bb700 = missing2nan(df.FLBBCD_BB_700_SCALED);

        # concatinate FLBBCD data into 1D arrays
        chla1d = cat(chla1d, chla, dims=1);
        cdom1d = cat(cdom1d, cdom, dims=1);
        bb1d = cat(bb1d, bb700, dims=1);

        # changing missing in AD2CP to NaN
        ad2cp_time = missing2nan(df.AD2CP_TIME)
        ad2cp_heading = missing2nan(df.AD2CP_HEADING)
        ad2cp_pitch = missing2nan(df.AD2CP_PITCH)
        ad2cp_roll = missing2nan(df.AD2CP_ROLL)
        ad2cp_pressure = missing2nan(df.AD2CP_PRESSURE)
        ad2cp_alt = missing2nan(df.AD2CP_ALT)
        ad2cp_v1_ch1 = missing2nan(df.AD2CP_V1_CH1)
        ad2cp_v2_ch1 = missing2nan(df.AD2CP_V2_CH1)
        ad2cp_v3_ch1 = missing2nan(df.AD2CP_V3_CH1)
        ad2cp_v4_ch1 = missing2nan(df.AD2CP_V4_CH1)
        ad2cp_v1_ch2 = missing2nan(df.AD2CP_V1_CH2)
        ad2cp_v2_ch2 = missing2nan(df.AD2CP_V2_CH2)
        ad2cp_v3_ch2 = missing2nan(df.AD2CP_V3_CH2)
        ad2cp_v4_ch2 = missing2nan(df.AD2CP_V4_CH2)
        ad2cp_v1_ch3 = missing2nan(df.AD2CP_V1_CH3)
        ad2cp_v2_ch3 = missing2nan(df.AD2CP_V2_CH3)
        ad2cp_v3_ch3 = missing2nan(df.AD2CP_V3_CH3)
        ad2cp_v4_ch3 = missing2nan(df.AD2CP_V4_CH3)
        ad2cp_v1_ch4 = missing2nan(df.AD2CP_V1_CH4)
        ad2cp_v2_ch4 = missing2nan(df.AD2CP_V2_CH4)
        ad2cp_v3_ch4 = missing2nan(df.AD2CP_V3_CH4)
        ad2cp_v4_ch4 = missing2nan(df.AD2CP_V4_CH4)
        ad2cp_v1_ch5 = missing2nan(df.AD2CP_V1_CH5)
        ad2cp_v2_ch5 = missing2nan(df.AD2CP_V2_CH5)
        ad2cp_v3_ch5 = missing2nan(df.AD2CP_V3_CH5)
        ad2cp_v4_ch5 = missing2nan(df.AD2CP_V4_CH5)
        ad2cp_v1_ch6 = missing2nan(df.AD2CP_V1_CH6)
        ad2cp_v2_ch6 = missing2nan(df.AD2CP_V2_CH6)
        ad2cp_v3_ch6 = missing2nan(df.AD2CP_V3_CH6)
        ad2cp_v4_ch6 = missing2nan(df.AD2CP_V4_CH6)
            
        # create profile data structure (organized by yo's) for NAV, CTD, and FLBBCD data
        push!(nav, NAV(t, z, lon, lat));
        push!(ctd, LEGATO(t, p, temp, cond, condtemp, salt));
        push!(flbbcd, FLBBCD(t, chla, cdom, bb700));
        push!(ad2cp, AD2CP(ad2cp_time, ad2cp_heading, ad2cp_pitch, ad2cp_roll, ad2cp_pressure, ad2cp_alt, ad2cp_v1_ch1, ad2cp_v2_ch1, ad2cp_v3_ch1, ad2cp_v4_ch1, ad2cp_v1_ch2, ad2cp_v2_ch2, ad2cp_v3_ch2, ad2cp_v4_ch2, ad2cp_v1_ch3, ad2cp_v2_ch3, ad2cp_v3_ch3, ad2cp_v4_ch3, ad2cp_v1_ch4, ad2cp_v2_ch4, ad2cp_v3_ch4, ad2cp_v4_ch4, ad2cp_v1_ch5, ad2cp_v2_ch5, ad2cp_v3_ch5, ad2cp_v4_ch5, ad2cp_v1_ch6, ad2cp_v2_ch6, ad2cp_v3_ch6, ad2cp_v4_ch6));
    end

    # storing 1d data into NAV, CTD, and FLBBCD data structures
    nav1d = NAV(t1d, z1d, lon1d, lat1d);
    ctd1d = LEGATO(t1d, p1d, temp1d, cond1d, condtemp1d, salt1d);
    flbbcd1d = FLBBCD(t1d, chla1d, cdom1d, bb1d);
    ad2cp1d = AD2CP(ad2cp_time1d, ad2cp_heading1d, ad2cp_pitch1d, ad2cp_roll1d, ad2cp_pressure1d, ad2cp_alt1d, ad2cp_v1_ch1_1d, ad2cp_v2_ch1_1d, ad2cp_v3_ch1_1d, ad2cp_v4_ch1_1d, ad2cp_v1_ch2_1d, ad2cp_v2_ch2_1d, ad2cp_v3_ch2_1d, ad2cp_v4_ch2_1d, ad2cp_v1_ch3_1d, ad2cp_v2_ch3_1d, ad2cp_v3_ch3_1d, ad2cp_v4_ch3_1d, ad2cp_v1_ch4_1d, ad2cp_v2_ch4_1d, ad2cp_v3_ch4_1d, ad2cp_v4_ch4_1d, ad2cp_v1_ch5_1d, ad2cp_v2_ch5_1d, ad2cp_v3_ch5_1d, ad2cp_v4_ch5_1d, ad2cp_v1_ch6_1d, ad2cp_v2_ch6_1d, ad2cp_v3_ch6_1d, ad2cp_v4_ch6_1d);

    # combinating NAV, CTD, and FLBBCD data into one glider data structure
    glider = SeaExplorer(nav, ctd, flbbcd, ad2cp, nav1d, ctd1d, flbbcd1d, ad2cp1d);

    return glider
end

# setting src and data directory paths
srcdir = "/Users/gong/GitHub/jlglider/";
dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";

# define dataset loading parameters
#project = "maracoos"
#deploydate = "20220311"
#suffix = "data"

project = "NORSE"
deploydate = "20220818"
suffix = "test"

glidername = "sea064"
mission = "34"

# define data load location
datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
navdir = datadir * "nav/";
scidir = datadir * "science/";

glider = load_SEAdata_raw(glidername, mission, navdir, scidir);
glider_rt = load_SEAdata_rt(glidername, mission, navdir, scidir);
