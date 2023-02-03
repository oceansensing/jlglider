# functions for working with AD2CP data
# gong@vims.edu 2022-10-24

module seaexplorer_functions

using Glob, DataFrames, CSV, Dates, Missings

import seaexplorer_types: NAV_RT, PLD_RT

function missing2nan(varin)
    varout = Float64.(collect(Missings.replace(varin, NaN)));
end

function cleanTime(varin)
    badtind = findall(varin .<= Dates.DateTime(2022,10,21,12,0,0));
    varout = varin;
    varout[badtind] .= NaN;
end

function cleanAD2CPtime(varin, varalt; systime = 1)

    if systime != 1
        varout = collect(copy(varin));
        badtind = findall(varin .== "00000000 00:00:00")
        if length(badtind) > 0
            for i = 1:length(badtind)
                tstr = varalt[badtind[i]];
                varout[badtind[i]] = tstr[4:5] * tstr[1:2] * tstr[9:10] * " " * tstr[12:19]
            end
        end
    elseif systime == 1
        varout = collect(varalt);
    end
    return varout
end

function cleanEPS(epsin)
    epsout = missing2nan(deepcopy(epsin));
    badind = findall((epsin .>= 1e-4) .|| (epsin .<= 1e-13));
    epsout[badind] .= NaN;
    epsout = convert(Vector{Float64}, epsout)
    return epsout;
end

function cleanTemp(varin)
    varout = missing2nan(deepcopy(varin));
    badind = findall((varin .>= 100.0) .|| (varin .<= -3));
    varout[badind] .= NaN;
    varout = convert(Vector{Float64}, varout);
    return varout;
end

function cleanSalt(varin)
    varout = missing2nan(deepcopy(varin));
    badind = findall((varin .>= 42.0) .|| (varin .<= 1.0));
    varout[badind] .= NaN;
    varout = convert(Vector{Float64}, varout);
    return varout;
end

function cleanPress(varin)
    varout = missing2nan(deepcopy(varin));
    badind = findall((varin .>= 7500.0) .|| (varin .< 0.0));
    varout[badind] .= NaN;
    varout = convert(Vector{Float64}, varout);
    return varout;
end

function clean9999(varin)
    varout = missing2nan(deepcopy(varin));
    badind = findall(varin .>= 7500.0);
    varout[badind] .= NaN;
    varout = convert(Vector{Float64}, varout);
    return varout;
end

function cleanAD2CP(varin)
    if typeof(collect(varin)) != Vector{Float64}
        varout = parse.(Float64, varin);
    else
        varout = convert.(Float64, varin);
    end
    #varout = deepcopy(varin);
    badind = findall(varout .<= -9000.0);
    varout[badind] .= NaN;
    #varout = convert(Vector{Float64}, varout);
    return varout;
end

#function cleanAD2CP(varin)
#end

function cleanFLBBCDchl(varin)
    varout = missing2nan(deepcopy(varin));
    badind = findall((varin .>= 50.0) .|| (varin .<= 0.0));
    varout[badind] .= NaN;
    varout = convert(Vector{Float64}, varout);
    return varout;
end

function cleanFLBBCDbb700(varin)
    varout = missing2nan(deepcopy(varin));
    badind = findall((varin .>= 1.0) .|| (varin .<= -1.0));
    varout[badind] .= NaN;
    varout = convert(Vector{Float64}, varout);
    return varout;
end

function cleanFLBBCDcdom(varin)
    varout = missing2nan(deepcopy(varin));
    badind = findall((varin .>= 50.0) .|| (varin .<= -2.0));
    varout[badind] .= NaN;
    varout = convert(Vector{Float64}, varout);
    return varout;
end

function load_NAV(gliderSN::Int, mission::Int, navdir::String, dataflag::Int)
    nav_rt = NAV_RT[];

    #missionroot = uppercase(glidername) * "." * mission;
    #gliroot_rt = missionroot * "." * "gli.sub.";
    #glilist_rt = Glob.glob(gliroot_rt * "*", navdir);

    glilist = Glob.glob( "*" * string(gliderSN; pad=3) * "." * string(mission) * ".gli.sub.*", navdir);

    # separating '.all' from '.###' files
    #glilist_suffix = [];
    yos = [];
    yolist = [];
    allindx = 0;
    for i = 1:length(glilist)
        suffix2 = glilist[i][end-1:end]
        fnlen = length(glilist[i]) - length(navdir);
        if suffix2 == "gz"
            if fnlen == 22  # sea064.37.gli.sub.1.gz
                yo = parse(Int,glilist[i][end-3:end-3]);
            elseif fnlen == 23 # sea064.37.gli.sub.10.gz
                yo = parse(Int,glilist[i][end-4:end-3]);
            elseif fnlen == 24 # sea064.37.gli.sub.100.gz
                yo = parse(Int,glilist[i][end-5:end-3]);
            elseif fnlen == 25 # sea064.37.gli.sub.1000.gz
                yo = parse(Int,glilist[i][end-6:end-3]);
            end
            yos = push!(yos, yo);
            yolist = push!(yolist, glilist[i]);
        elseif (suffix2 != "gz" && suffix2 != "ll")
            if fnlen == 19  # sea064.37.gli.sub.1
                yo = parse(Int,glilist[i][end:end]);
            elseif fnlen == 20 # sea064.37.gli.sub.10
                yo = parse(Int,glilist[i][end-1:end]);
            elseif fnlen == 21 # sea064.37.gli.sub.100
                yo = parse(Int,glilist[i][end-2:end]);
            elseif fnlen == 22 # sea064.37.gli.sub.1000
                yo = parse(Int,glilist[i][end-3:end]);
            end
            yos = push!(yos, yo);
            yolist = push!(yolist, glilist[i]);
        elseif suffix2 == "ll"
            if fnlen == 21 # SEA064.37.gli.sub.all
                allindx = i;
            end
        end
        #glilist_suffix = push!(glilist_suffix, glilist[i][end-2:end]);
    end
    #yolist = findall(glilist_suffix .!= "all");
    #allindx = findall(glilist_suffix .== "all");
    yosi = sortperm(Int.(yos));
    yos = yos[yosi];
    yolist = yolist[yosi];

    if dataflag == 1
        glilist = [glilist[allindx]];
    else
        glilist = yolist;
    end

    #Glob.glob("SEA064.37.pld1.sub.*", navdir);

    # time data format
    timeformat = "dd/mm/yyyy HH:MM:SS.sss"
    
    # initiate NAV_RT 1-D variables
    yo1d = [];
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
    for i = 1:length(glilist)
        #yostring = glilist[i][end-2:end];

        if dataflag == 1
        #    yo = parse.(Int, glilist_suffix[yolist]);
            yo = yos;
        else
        #    yo = [parse(Int, yostring)];
            yo = [yos[i]];
        end

        # load nav files, handle .gz if they are compressed
        ##navfilepath = navdir * gliroot_rt * string(i, pad=3); 
        #navfilepath = navdir * gliroot_rt * yostring * ".gz"; 
        #if isfile(navfilepath) != true
        #    navfilepath = navdir * gliroot_rt * yostring; 
        #end

        navfilepath = glilist[i];
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

        yo1d = cat(yo1d, yo, dims=1);
        t1d = cat(t1d, t, dims=1);
        z1d = cat(z1d, z, dims=1);
        lon1d = cat(lon1d, lon, dims=1);
        lat1d = cat(lat1d, lat, dims=1);
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

        push!(nav_rt, NAV_RT(yo, t, z, lon, lat, NavState, SecurityLevel, Heading, Declination, Pitch, Roll, DeadReckoning, DesiredH, BallastCmd, BallastPos, LinCmd, LinPos, AngCmd, AngPos, Voltage, Altitude));
    end #for
    nav1d_rt = NAV_RT(yo1d, t1d, z1d, lon1d, lat1d, NavState1d, SecurityLevel1d, Heading1d, Declination1d, Pitch1d, Roll1d, DeadReckoning1d, DesiredH1d, BallastCmd1d, BallastPos1d, LinCmd1d, LinPos1d, AngCmd1d, AngPos1d, Voltage1d, Altitude1d);

    # combinating NAV_RT and PLD_RT data into one glider data structure
    #gliderRT = SeaExplorerRT(nav_rt, pld_rt, nav1d_rt, pld1d_rt);

    return nav_rt, nav1d_rt
end

function load_PLD(gliderSN::Int, mission::Int, scidir::String, dataflag::Int)
    pld_rt = PLD_RT[];

    if dataflag < 2
        datatype = "sub"
    else
        datatype = "raw"
    end

    #missionroot = uppercase(glidername) * "." * mission;
    #pldroot_rt = missionroot * "." * "pld1.sub.";
    #pldlist_rt = Glob.glob(pldroot_rt * "*", scidir);

    pldlist = Glob.glob("*" * string(gliderSN; pad=3) * "." * string(mission) * ".pld1." * datatype * ".*", scidir);

    #pldlist_suffix =[];
    #for i = 1:length(pldlist_rt)
    #    pldlist_suffix = push!(pldlist_suffix, pldlist_rt[i][end-2:end]);
    #end
    #yolist = findall(pldlist_suffix .!= "all");
    #allindx = findall(pldlist_suffix .== "all");

    yos = [];
    yolist = [];
    allindx = 0;
    for i = 1:length(pldlist)
        suffix2 = pldlist[i][end-1:end]
        fnlen = length(pldlist[i]) - length(scidir);
        if suffix2 == "gz"
            if fnlen == 22+1  # sea064.37.pld1.sub.1.gz
                yo = parse(Int,pldlist[i][end-3:end-3]);
            elseif fnlen == 23+1 # sea064.37.pld1.sub.10.gz
                yo = parse(Int,pldlist[i][end-4:end-3]);
            elseif fnlen == 24+1 # sea064.37.pld1.sub.100.gz
                yo = parse(Int,pldlist[i][end-5:end-3]);
            elseif fnlen == 25+1 # sea064.37.pld1.sub.1000.gz
                yo = parse(Int,pldlist[i][end-6:end-3]);
            end
            yos = push!(yos, yo);
            yolist = push!(yolist, pldlist[i]);
        elseif (suffix2 != "gz" && suffix2 != "ll")
            if fnlen == 19+1  # sea064.37.pld1.sub.1
                yo = parse(Int,pldlist[i][end:end]);
            elseif fnlen == 20+1 # sea064.37.pld1.sub.10
                yo = parse(Int,pldlist[i][end-1:end]);
            elseif fnlen == 21+1 # sea064.37.pld1.sub.100
                yo = parse(Int,pldlist[i][end-2:end]);
            elseif fnlen == 22+1 # sea064.37.pld1.sub.1000
                yo = parse(Int,pldlist[i][end-3:end]);
            end
            yos = push!(yos, yo);
            yolist = push!(yolist, pldlist[i]);
        elseif suffix2 == "ll"
            if fnlen == 21+1 # SEA064.37.pld1.sub.all
                allindx = i;
            end
        end
        #glilist_suffix = push!(glilist_suffix, glilist[i][end-2:end]);
    end
    #yolist = findall(glilist_suffix .!= "all");
    #allindx = findall(glilist_suffix .== "all");
    yosi = sortperm(Int.(yos));
    yos = yos[yosi];
    yolist = yolist[yosi];

    if dataflag == 1
        pldlist = [pldlist[allindx]];
    else
        pldlist = yolist;
    end

    #if dataflag == 1
    #    pldlist_rt = pldlist_rt[allindx];
    #else
    #    pldlist_rt = pldlist_rt[yolist];
    #end

    # time data format
    timeformat = "dd/mm/yyyy HH:MM:SS.sss"
    ad2cptimeformat = "mmddyy HH:MM:SS"
    
    # initiate PLD_RT 1-D variables
    yo1d = [];
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
    ad2cp_v1_cn1_1d = [];
    ad2cp_v2_cn1_1d = [];
    ad2cp_v3_cn1_1d = [];
    ad2cp_v4_cn1_1d = [];
    ad2cp_v1_cn2_1d = [];
    ad2cp_v2_cn2_1d = [];
    ad2cp_v3_cn2_1d = [];
    ad2cp_v4_cn2_1d = [];
    ad2cp_v1_cn3_1d = [];
    ad2cp_v2_cn3_1d = [];
    ad2cp_v3_cn3_1d = [];
    ad2cp_v4_cn3_1d = [];
    ad2cp_v1_cn4_1d = [];
    ad2cp_v2_cn4_1d = [];
    ad2cp_v3_cn4_1d = [];
    ad2cp_v4_cn4_1d = [];
    ad2cp_v1_cn5_1d = [];
    ad2cp_v2_cn5_1d = [];
    ad2cp_v3_cn5_1d = [];
    ad2cp_v4_cn5_1d = [];
    ad2cp_v1_cn6_1d = [];
    ad2cp_v2_cn6_1d = [];
    ad2cp_v3_cn6_1d = [];
    ad2cp_v4_cn6_1d = [];
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
    ad2cp_Unorth_1d = [];
    ad2cp_Ueast_1d = [];
    ad2cp_Utot_1d = [];
    ad2cp_Udir_1d = [];
    ad2cp_qf_1d = [];

    # loop through the list of realtime payload (pld) data files
    for i = 1:length(pldlist)
        #yostring = pldlist[i][end-2:end];

        # separating '.all' from '.###' files
        if dataflag == 1
        #    yo = parse.(Int, pldlist_suffix[yolist]);
            yo = yos
        else
        #    yo = [parse(Int, yostring)];
            yo = [yos[i]];
        end

        # load science files, handle .gz if they are compressed
        ##pldfilepath = scidir * pldroot_rt * string(i, pad=3); 
        #pldfilepath = scidir * pldroot_rt * yostring * ".gz";
        #if isfile(pldfilepath) != true
        #    pldfilepath = scidir * pldroot_rt * yostring;
        #end

        pldfilepath = pldlist[i];
        print(pldfilepath * "\n")
        df = CSV.read(pldfilepath, header=1, delim=";", DataFrame, buffer_in_memory=true);

        # extract location data from data frame
        t = DateTime.(df.PLD_REALTIMECLOCK, timeformat);
        z = missing2nan(df.NAV_DEPTH);
        navlon = df.NAV_LONGITUDE;
        navlat = df.NAV_LATITUDE;
        lon = trunc.(navlon ./ 100) + (navlon .- trunc.(navlon ./ 100)*100) / 60;
        lat = trunc.(navlat ./ 100) + (navlat .- trunc.(navlat ./ 100)*100) / 60;
        lon = missing2nan(lon);
        lat = missing2nan(lat);
        nav_resource = df.NAV_RESOURCE;
        ad2cp_time = DateTime.(cleanAD2CPtime(collect(df.AD2CP_TIME), collect(df.PLD_REALTIMECLOCK)), ad2cptimeformat) .+ Dates.Year(2000);
        ad2cp_heading = df.AD2CP_HEADING;
        ad2cp_pitch = df.AD2CP_PITCH;
        ad2cp_roll = df.AD2CP_ROLL;
        ad2cp_pressure = df.AD2CP_PRESSURE;
        ad2cp_alt = df.AD2CP_ALT;
        ad2cp_v1_cn1 = cleanAD2CP(df.AD2CP_V1_CN1);
        ad2cp_v2_cn1 = cleanAD2CP(df.AD2CP_V2_CN1);
        ad2cp_v3_cn1 = cleanAD2CP(df.AD2CP_V3_CN1);
        ad2cp_v4_cn1 = cleanAD2CP(df.AD2CP_V4_CN1);
        ad2cp_v1_cn2 = cleanAD2CP(df.AD2CP_V1_CN2);
        ad2cp_v2_cn2 = cleanAD2CP(df.AD2CP_V2_CN2);
        ad2cp_v3_cn2 = cleanAD2CP(df.AD2CP_V3_CN2);
        ad2cp_v4_cn2 = cleanAD2CP(df.AD2CP_V4_CN2);
        ad2cp_v1_cn3 = cleanAD2CP(df.AD2CP_V1_CN3);
        ad2cp_v2_cn3 = cleanAD2CP(df.AD2CP_V2_CN3);
        ad2cp_v3_cn3 = cleanAD2CP(df.AD2CP_V3_CN3);
        ad2cp_v4_cn3 = cleanAD2CP(df.AD2CP_V4_CN3);
        ad2cp_v1_cn4 = cleanAD2CP(df.AD2CP_V1_CN4);
        ad2cp_v2_cn4 = cleanAD2CP(df.AD2CP_V2_CN4);
        ad2cp_v3_cn4 = cleanAD2CP(df.AD2CP_V3_CN4);
        ad2cp_v4_cn4 = cleanAD2CP(df.AD2CP_V4_CN4);
        ad2cp_v1_cn5 = cleanAD2CP(df.AD2CP_V1_CN5);
        ad2cp_v2_cn5 = cleanAD2CP(df.AD2CP_V2_CN5);
        ad2cp_v3_cn5 = cleanAD2CP(df.AD2CP_V3_CN5);
        ad2cp_v4_cn5 = cleanAD2CP(df.AD2CP_V4_CN5);
        ad2cp_v1_cn6 = cleanAD2CP(df.AD2CP_V1_CN6);
        ad2cp_v2_cn6 = cleanAD2CP(df.AD2CP_V2_CN6);
        ad2cp_v3_cn6 = cleanAD2CP(df.AD2CP_V3_CN6);
        ad2cp_v4_cn6 = cleanAD2CP(df.AD2CP_V4_CN6);
        flbbcd_chl_count = df.FLBBCD_CHL_COUNT;
        flbbcd_chl_scaled = df.FLBBCD_CHL_SCALED;
        flbbcd_bb_700_count = df.FLBBCD_BB_700_COUNT;
        flbbcd_bb_700_scaled = df.FLBBCD_BB_700_SCALED;
        flbbcd_cdom_count = df.FLBBCD_CDOM_COUNT;
        flbbcd_cdom_scaled = df.FLBBCD_CDOM_SCALED;
        legato_conductivity = df.LEGATO_CONDUCTIVITY;
        legato_temperature = df.LEGATO_TEMPERATURE;
        legato_pressure = df.LEGATO_PRESSURE;
        legato_salinity = df.LEGATO_SALINITY;
        legato_condtemp = df.LEGATO_CONDTEMP;
        mr1000g_t1_avg = df."MR1000G-RDL_T1_AVG";
        mr1000g_t2_avg = df."MR1000G-RDL_T2_AVG";
        mr1000g_sh1_std = df."MR1000G-RDL_SH1_STD";
        mr1000g_sh2_std = df."MR1000G-RDL_SH2_STD";
        mr1000g_press_avg = df."MR1000G-RDL_PRESS_AVG";
        mr1000g_incly_avg = df."MR1000G-RDL_INCLY_AVG";
        mr1000g_eps1 = df."MR1000G-RDL_EPS1";
        mr1000g_qc1 = df."MR1000G-RDL_QC1";
        mr1000g_eps2 = df."MR1000G-RDL_EPS2";
        mr1000g_qc2 = df."MR1000G-RDL_QC2";

        tmpvar = Array{Float64,1}(undef, length(ad2cp_alt));
        tmpvar .= NaN;

        if dataflag == 1
            ad2cp_Unorth = df.AD2CP_Unorth_c;
            ad2cp_Ueast = df.AD2CP_Ueast_c;
            ad2cp_Utot = df.AD2CP_Utot_c;
            ad2cp_Udir = df.AD2CP_Udir_c;
            ad2cp_qf = df.AD2CP_QF_c;
        else
            ad2cp_Unorth = tmpvar;
            ad2cp_Ueast = tmpvar;
            ad2cp_Utot = tmpvar;
            ad2cp_Udir = tmpvar;
            ad2cp_qf = tmpvar;
        end

        yo1d = cat(yo1d, yo[1], dims = 1);
        t1d = cat(t1d, t, dims = 1);
        z1d = cat(z1d, z, dims = 1);
        lon1d = cat(lon1d, lon, dims = 1);
        lat1d = cat(lat1d, lat, dims = 1);
        nav_resource1d = cat(nav_resource1d, nav_resource, dims = 1);
        ad2cp_time1d = cat(ad2cp_time1d, ad2cp_time, dims = 1);
        ad2cp_heading1d = cat(ad2cp_heading1d, ad2cp_heading, dims = 1);
        ad2cp_pitch1d = cat(ad2cp_pitch1d, ad2cp_pitch, dims = 1);
        ad2cp_roll1d = cat(ad2cp_roll1d, ad2cp_roll, dims = 1);
        ad2cp_pressure1d = cat(ad2cp_pressure1d, ad2cp_pressure, dims = 1);
        ad2cp_alt1d = cat(ad2cp_alt1d, ad2cp_alt, dims = 1);
        ad2cp_v1_cn1_1d = cat(ad2cp_v1_cn1_1d, ad2cp_v1_cn1, dims = 1);
        ad2cp_v2_cn1_1d = cat(ad2cp_v2_cn1_1d, ad2cp_v2_cn1, dims = 1);
        ad2cp_v3_cn1_1d = cat(ad2cp_v3_cn1_1d, ad2cp_v3_cn1, dims = 1);
        ad2cp_v4_cn1_1d = cat(ad2cp_v4_cn1_1d, ad2cp_v4_cn1, dims = 1);
        ad2cp_v1_cn2_1d = cat(ad2cp_v1_cn2_1d, ad2cp_v1_cn2, dims = 1);
        ad2cp_v2_cn2_1d = cat(ad2cp_v2_cn2_1d, ad2cp_v2_cn2, dims = 1);
        ad2cp_v3_cn2_1d = cat(ad2cp_v3_cn2_1d, ad2cp_v3_cn2, dims = 1);
        ad2cp_v4_cn2_1d = cat(ad2cp_v4_cn2_1d, ad2cp_v4_cn2, dims = 1);
        ad2cp_v1_cn3_1d = cat(ad2cp_v1_cn3_1d, ad2cp_v1_cn3, dims = 1);
        ad2cp_v2_cn3_1d = cat(ad2cp_v2_cn3_1d, ad2cp_v2_cn3, dims = 1);
        ad2cp_v3_cn3_1d = cat(ad2cp_v3_cn3_1d, ad2cp_v3_cn3, dims = 1);
        ad2cp_v4_cn3_1d = cat(ad2cp_v4_cn3_1d, ad2cp_v4_cn3, dims = 1);
        ad2cp_v1_cn4_1d = cat(ad2cp_v1_cn4_1d, ad2cp_v1_cn4, dims = 1);
        ad2cp_v2_cn4_1d = cat(ad2cp_v2_cn4_1d, ad2cp_v2_cn4, dims = 1);
        ad2cp_v3_cn4_1d = cat(ad2cp_v3_cn4_1d, ad2cp_v3_cn4, dims = 1);
        ad2cp_v4_cn4_1d = cat(ad2cp_v4_cn4_1d, ad2cp_v4_cn4, dims = 1);
        ad2cp_v1_cn5_1d = cat(ad2cp_v1_cn5_1d, ad2cp_v1_cn5, dims = 1);
        ad2cp_v2_cn5_1d = cat(ad2cp_v2_cn5_1d, ad2cp_v2_cn5, dims = 1);
        ad2cp_v3_cn5_1d = cat(ad2cp_v3_cn5_1d, ad2cp_v3_cn5, dims = 1);
        ad2cp_v4_cn5_1d = cat(ad2cp_v4_cn5_1d, ad2cp_v4_cn5, dims = 1);
        ad2cp_v1_cn6_1d = cat(ad2cp_v1_cn6_1d, ad2cp_v1_cn6, dims = 1);
        ad2cp_v2_cn6_1d = cat(ad2cp_v2_cn6_1d, ad2cp_v2_cn6, dims = 1);
        ad2cp_v3_cn6_1d = cat(ad2cp_v3_cn6_1d, ad2cp_v3_cn6, dims = 1);
        ad2cp_v4_cn6_1d = cat(ad2cp_v4_cn6_1d, ad2cp_v4_cn6, dims = 1);
        flbbcd_chl_count_1d = cat(flbbcd_chl_count_1d, flbbcd_chl_count, dims = 1);
        flbbcd_chl_scaled_1d = cat(flbbcd_chl_scaled_1d, flbbcd_chl_scaled, dims = 1);
        flbbcd_bb_700_count_1d = cat(flbbcd_bb_700_count_1d, flbbcd_bb_700_count, dims = 1);
        flbbcd_bb_700_scaled_1d = cat(flbbcd_bb_700_scaled_1d, flbbcd_bb_700_scaled, dims = 1);
        flbbcd_cdom_count_1d = cat(flbbcd_cdom_count_1d, flbbcd_cdom_count, dims = 1);
        flbbcd_cdom_scaled_1d = cat(flbbcd_cdom_scaled_1d, flbbcd_cdom_scaled, dims = 1);
        legato_conductivity_1d = cat(legato_conductivity_1d, legato_conductivity, dims = 1);
        legato_temperature_1d = cat(legato_temperature_1d, legato_temperature, dims = 1);
        legato_pressure_1d = cat(legato_pressure_1d, legato_pressure, dims = 1);
        legato_salinity_1d = cat(legato_salinity_1d, legato_salinity, dims = 1);
        legato_condtemp_1d = cat(legato_condtemp_1d, legato_condtemp, dims = 1);
        mr1000g_t1_avg_1d = cat(mr1000g_t1_avg_1d, mr1000g_t1_avg, dims = 1);
        mr1000g_t2_avg_1d = cat(mr1000g_t2_avg_1d, mr1000g_t2_avg, dims = 1);
        mr1000g_sh1_std_1d = cat(mr1000g_sh1_std_1d, mr1000g_sh1_std, dims = 1);
        mr1000g_sh2_std_1d = cat(mr1000g_sh2_std_1d, mr1000g_sh2_std, dims = 1);
        mr1000g_press_avg_1d = cat(mr1000g_press_avg_1d, mr1000g_press_avg, dims = 1);
        mr1000g_incly_avg_1d = cat(mr1000g_incly_avg_1d, mr1000g_incly_avg, dims = 1);
        mr1000g_eps1_1d = cat(mr1000g_eps1_1d, mr1000g_eps1, dims = 1);
        mr1000g_qc1_1d = cat(mr1000g_qc1_1d, mr1000g_qc1, dims = 1);
        mr1000g_eps2_1d = cat(mr1000g_eps2_1d, mr1000g_eps2, dims = 1);
        mr1000g_qc2_1d = cat(mr1000g_qc2_1d, mr1000g_qc2, dims = 1);

        ad2cp_Unorth_1d = cat(ad2cp_Unorth_1d, ad2cp_Unorth, dims = 1);
        ad2cp_Ueast_1d = cat(ad2cp_Ueast_1d, ad2cp_Ueast, dims = 1);
        ad2cp_Utot_1d = cat(ad2cp_Utot_1d, ad2cp_Utot, dims = 1);
        ad2cp_Udir_1d = cat(ad2cp_Udir_1d, ad2cp_Udir, dims = 1);
        ad2cp_qf_1d = cat(ad2cp_qf_1d, ad2cp_qf, dims = 1);
        
        push!(pld_rt, PLD_RT(yo, t, z, lon, lat, nav_resource, ad2cp_time, ad2cp_heading, ad2cp_pitch, ad2cp_roll, ad2cp_pressure, ad2cp_alt,  ad2cp_v1_cn1, ad2cp_v2_cn1, ad2cp_v3_cn1, ad2cp_v4_cn1, ad2cp_v1_cn2, ad2cp_v2_cn2, ad2cp_v3_cn2, ad2cp_v4_cn2, ad2cp_v1_cn3, ad2cp_v2_cn3, ad2cp_v3_cn3, ad2cp_v4_cn3, ad2cp_v1_cn4, ad2cp_v2_cn4, ad2cp_v3_cn4, ad2cp_v4_cn4, ad2cp_v1_cn5, ad2cp_v2_cn5, ad2cp_v3_cn5, ad2cp_v4_cn5, ad2cp_v1_cn6, ad2cp_v2_cn6, ad2cp_v3_cn6, ad2cp_v4_cn6, flbbcd_chl_count, flbbcd_chl_scaled, flbbcd_bb_700_count, flbbcd_bb_700_scaled, flbbcd_cdom_count, flbbcd_cdom_scaled, legato_conductivity, legato_temperature, legato_pressure, legato_salinity, legato_condtemp, mr1000g_t1_avg, mr1000g_t2_avg, mr1000g_sh1_std, mr1000g_sh2_std, mr1000g_press_avg, mr1000g_incly_avg, mr1000g_eps1, mr1000g_qc1, mr1000g_eps2, mr1000g_qc2, ad2cp_Unorth, ad2cp_Ueast, ad2cp_Utot, ad2cp_Udir, ad2cp_qf));    
    end #for
    pld1d_rt = PLD_RT(yo1d, t1d, z1d, lon1d, lat1d, nav_resource1d, ad2cp_time1d, ad2cp_heading1d, ad2cp_pitch1d, ad2cp_roll1d, ad2cp_pressure1d, ad2cp_alt1d, ad2cp_v1_cn1_1d, ad2cp_v2_cn1_1d, ad2cp_v3_cn1_1d, ad2cp_v4_cn1_1d, ad2cp_v1_cn2_1d, ad2cp_v2_cn2_1d, ad2cp_v3_cn2_1d, ad2cp_v4_cn2_1d, ad2cp_v1_cn3_1d, ad2cp_v2_cn3_1d, ad2cp_v3_cn3_1d, ad2cp_v4_cn3_1d, ad2cp_v1_cn4_1d, ad2cp_v2_cn4_1d, ad2cp_v3_cn4_1d, ad2cp_v4_cn4_1d, ad2cp_v1_cn5_1d, ad2cp_v2_cn5_1d, ad2cp_v3_cn5_1d, ad2cp_v4_cn5_1d, ad2cp_v1_cn6_1d, ad2cp_v2_cn6_1d, ad2cp_v3_cn6_1d, ad2cp_v4_cn6_1d, flbbcd_chl_count_1d, flbbcd_chl_scaled_1d, flbbcd_bb_700_count_1d, flbbcd_bb_700_scaled_1d, flbbcd_cdom_count_1d, flbbcd_cdom_scaled_1d, legato_conductivity_1d, legato_temperature_1d, legato_pressure_1d, legato_salinity_1d, legato_condtemp_1d, mr1000g_t1_avg_1d, mr1000g_t2_avg_1d, mr1000g_sh1_std_1d, mr1000g_sh2_std_1d, mr1000g_press_avg_1d, mr1000g_incly_avg_1d, mr1000g_eps1_1d, mr1000g_qc1_1d, mr1000g_eps2_1d, mr1000g_qc2_1d, ad2cp_Unorth_1d, ad2cp_Ueast_1d, ad2cp_Utot_1d, ad2cp_Udir_1d, ad2cp_qf_1d);

    # combinating NAV_RT and PLD_RT data into one glider data structure
    #gliderRT = SeaExplorerRT(nav_rt, pld_rt, nav1d_rt, pld1d_rt);

    return pld_rt, pld1d_rt
end

end