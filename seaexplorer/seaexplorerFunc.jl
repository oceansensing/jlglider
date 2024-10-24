# functions for working with AD2CP data
# gong@vims.edu 2022-10-24

module seaexplorerFunc

include("seaexplorerType.jl");
include("gsw_c2po.jl");

using Glob, DataFrames, CSV, Dates, Missings, NaNMath, Interpolations, YAML
using NCDatasets
using GibbsSeaWater, MAT
import .seaexplorerType: NAV_RT, PLD_RT, SeaExplorerData
import .gsw_c2po: sigma0_from_t_sp, spice0_from_t_sp, N2_from_t_sp

# https://discourse.julialang.org/t/indices-of-intersection-of-two-arrays/23043/20
function intersectalajulia2(a,b)
    ia = findall(in(b), a)
    ib = findall(in(view(a,ia)), b)
    return unique(view(a,ia)), ia, ib[indexin(view(a,ia), view(b,ib))]
end

function yearday(xdt::DateTime)
    yday = Dates.dayofyear(xdt);
    seconds_in_day = 86400;
    ydayfrac = yday + (Dates.hour(xdt) * 3600 .+ Dates.minute(xdt) * 60 .+ Dates.second(xdt)) ./ seconds_in_day;
    return ydayfrac;
end

function yearday(unix_t::AbstractFloat)
    xdt = unix2datetime(unix_t); 
    yday = Dates.dayofyear(xdt);
    seconds_in_day = 86400;
    ydayfrac = yday + (Dates.hour(xdt) * 3600 + Dates.minute(xdt) * 60 .+ Dates.second(xdt)) ./ seconds_in_day;
    return ydayfrac;
end

function datetick(unix_t)
    x = unix_t
    xdt = unix2datetime.(unix_t); 
    yday = Dates.dayofyear.(xdt);
    uyday = unique(yday);
    hour = Dates.Hour.(xdt);
    minute = Dates.Minute.(xdt);
    #df = DateFormat("y-m-d");

    #tickind = Vector{Int64}(undef, length(uyday));
    tickind = [];
    for i = 1:length(uyday)
        t0 = findall((yday .== uyday[i]) .& (hour .== Hour(0)));
        if isempty(t0) == false
            push!(tickind, t0[1]);
        end
    end
    xtick = x[tickind];
    #xticklabel = string.(Dates.Date.(xdt[tickind]));
    xticklabel = [x[6:10] for x in string.(xdt[tickind])];
    return xdt, xtick, xticklabel    
end

function missing2nan(varin)
    varin = collect(varin);
    if (typeof(varin) == Vector{Union{Missing, Int64}}) | (typeof(varin) == Matrix{Union{Missing, Int64}})
        varout = Array{Float64}(undef,size(collect(varin)));
        varintypes = typeof.(varin);
        notmissind = findall(varintypes .!= Missing);
        missind = findall(varintypes .== Missing); 
        if isempty(notmissind) != true  
            varout[notmissind] .= Float64.(varin[notmissind]);
        end
        if isempty(missind) != true
            varout[missind] .= NaN;
        end
    elseif (typeof(varin) == Vector{Union{Missing, Float64}}) | (typeof(varin) == Matrix{Union{Missing, Float64}})
        varout = Float64.(collect(Missings.replace(varin, NaN)));
    elseif (typeof(varin) == Vector{Missing}) | (typeof(varin) == Matrix{Missing})
        varout = Array{Float64}(undef,size(collect(varin)));
        varout .= NaN; 
    else
        varout = varin;
    end

    return varout
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
    badind = findall((varin .>= 42.0) .|| (varin .<= 5.0));
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
    varin = missing2nan(varin);
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

function load_NAV(gliderSN::Int, missionID::Int, project::String, navdir::String, dataflag::String)
    nav_rt = NAV_RT[];

    #missionroot = uppercase(glidername) * "." * mission;
    #gliroot_rt = missionroot * "." * "gli.sub.";
    #glilist_rt = Glob.glob(gliroot_rt * "*", navdir);

    glilist = Glob.glob( "*" * string(gliderSN; pad=3) * "." * string(missionID) * ".gli.sub.*", navdir);

    # separating '.all' from '.###' files
    #glilist_suffix = [];
    yos = [];
    yolist = [];
    allindx = 1;
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
        elseif (suffix2 != "gz") & (suffix2 != "ll") & (suffix2 != "sv")
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
            yos = [-1];
        elseif suffix2 == "sv"
            if fnlen == 25 # SEA064.37.gli.sub.all.csv
                allindx = i;
            end
            yos = [-1];
        end
        #glilist_suffix = push!(glilist_suffix, glilist[i][end-2:end]);
    end
    #yolist = findall(glilist_suffix .!= "all");
    #allindx = findall(glilist_suffix .== "all");
    yosi = sortperm(Int.(yos));
    yos = yos[yosi];

    if !isempty(yolist)
        yolist = yolist[yosi];
    end

    if dataflag == "all"
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

        # load nav files, handle .gz if they are compressed
        ##navfilepath = navdir * gliroot_rt * string(i, pad=3); 
        #navfilepath = navdir * gliroot_rt * yostring * ".gz"; 
        #if isfile(navfilepath) != true
        #    navfilepath = navdir * gliroot_rt * yostring; 
        #end

        navfilepath = glilist[i];
        print(navfilepath * "\n")
        df = CSV.read(navfilepath, header=1, delim=";", DataFrame, buffer_in_memory=true);

        if dataflag == "all"
        #    yo = parse.(Int, glilist_suffix[yolist]);
            yo = df.YO_NUMBER;
        else
        #    yo = [parse(Int, yostring)];
            yo = [yos[i]];
        end
    
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

        push!(nav_rt, NAV_RT(gliderSN, missionID, project, yo, t, z, lon, lat, NavState, SecurityLevel, Heading, Declination, Pitch, Roll, DeadReckoning, DesiredH, BallastCmd, BallastPos, LinCmd, LinPos, AngCmd, AngPos, Voltage, Altitude));
    end #for
    nav1d_rt = NAV_RT(gliderSN, missionID, project, yo1d, t1d, z1d, lon1d, lat1d, NavState1d, SecurityLevel1d, Heading1d, Declination1d, Pitch1d, Roll1d, DeadReckoning1d, DesiredH1d, BallastCmd1d, BallastPos1d, LinCmd1d, LinPos1d, AngCmd1d, AngPos1d, Voltage1d, Altitude1d);

    # combinating NAV_RT and PLD_RT data into one glider data structure
    #gliderRT = SeaExplorerRT(nav_rt, pld_rt, nav1d_rt, pld1d_rt);

    return nav_rt, nav1d_rt
end

function load_PLD(gliderSN::Int, missionID::Int, project::String, scidir::String, dataflag::String)
    pld_rt = PLD_RT[];

    #if dataflag < 2
    #    datatype = "sub"
    #else
    #    datatype = "raw"
    #end

    if (dataflag == "sub") | (dataflag == "realtime") | (dataflag == "all")
        datatype = "sub"
    elseif (dataflag == "raw") | (dataflag == "delayed")
        datatype = "raw"
    end

    #missionroot = uppercase(glidername) * "." * mission;
    #pldroot_rt = missionroot * "." * "pld1.sub.";
    #pldlist_rt = Glob.glob(pldroot_rt * "*", scidir);

    pldlist = Glob.glob("*" * string(gliderSN; pad=3) * "." * string(missionID) * ".pld1." * datatype * ".*", scidir);

    #pldlist_suffix =[];
    #for i = 1:length(pldlist_rt)
    #    pldlist_suffix = push!(pldlist_suffix, pldlist_rt[i][end-2:end]);
    #end
    #yolist = findall(pldlist_suffix .!= "all");
    #allindx = findall(pldlist_suffix .== "all");

    yos = [];
    yolist = [];
    allindx = 1;
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
        elseif (suffix2 != "gz") & (suffix2 != "ll") & (suffix2 != "sv")
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
            yos = [-1];
        elseif suffix2 == "sv"
            if fnlen == 25 # SEA064.37.gli.sub.all.csv
                allindx = i;
            end
            yos = [-1];
        end
        #glilist_suffix = push!(glilist_suffix, glilist[i][end-2:end]);
    end
    #yolist = findall(glilist_suffix .!= "all");
    #allindx = findall(glilist_suffix .== "all");
    yosi = sortperm(Int.(yos));
    yos = yos[yosi];

    if !isempty(yolist)
        yolist = yolist[yosi];
    end

    if dataflag == "all"
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

        # load science files, handle .gz if they are compressed
        ##pldfilepath = scidir * pldroot_rt * string(i, pad=3); 
        #pldfilepath = scidir * pldroot_rt * yostring * ".gz";
        #if isfile(pldfilepath) != true
        #    pldfilepath = scidir * pldroot_rt * yostring;
        #end

        pldfilepath = pldlist[i];
        print(pldfilepath * "\n")
        df = CSV.read(pldfilepath, header=1, delim=";", DataFrame, buffer_in_memory=true);

        # separating '.all' from '.###' files
        if dataflag == "all"
        #    yo = parse.(Int, pldlist_suffix[yolist]);
            yo = df.YO_NUMBER;
        else
        #    yo = [parse(Int, yostring)];
            yo = [yos[i]];
        end
    
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
        #ad2cp_time = DateTime.(cleanAD2CPtime(collect(df.AD2CP_TIME), collect(df.PLD_REALTIMECLOCK)), ad2cptimeformat) .+ Dates.Year(2000);
        ad2cp_time = t;
        ad2cp_heading = missing2nan(df.AD2CP_HEADING);
        ad2cp_pitch = missing2nan(df.AD2CP_PITCH);
        ad2cp_roll = missing2nan(df.AD2CP_ROLL);
        ad2cp_pressure = missing2nan(df.AD2CP_PRESSURE);
        ad2cp_alt = missing2nan(df.AD2CP_ALT);
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

        #flbbcd_chl_count = missing2nan(df.FLBBCD_CHL_COUNT);
        #flbbcd_chl_scaled = missing2nan(df.FLBBCD_CHL_SCALED);
        #flbbcd_bb_700_count = missing2nan(df.FLBBCD_BB_700_COUNT);
        #flbbcd_bb_700_scaled = missing2nan(df.FLBBCD_BB_700_SCALED);
        #flbbcd_cdom_count = missing2nan(df.FLBBCD_CDOM_COUNT);
        #flbbcd_cdom_scaled = missing2nan(df.FLBBCD_CDOM_SCALED);

        flbbcd_chl_count = [];
        flbbcd_chl_scaled = [];
        flbbcd_bb_700_count = [];
        flbbcd_bb_700_scaled = [];
        flbbcd_cdom_count = [];
        flbbcd_cdom_scaled = [];

        legato_conductivity = missing2nan(df.LEGATO_CONDUCTIVITY);
        legato_temperature = missing2nan(df.LEGATO_TEMPERATURE);
        legato_pressure = missing2nan(df.LEGATO_PRESSURE);
        legato_salinity = missing2nan(df.LEGATO_SALINITY);
        legato_condtemp = missing2nan(df.LEGATO_CONDTEMP);
        mr1000g_t1_avg = missing2nan(df."MR1000G-RDL_T1_AVG");
        mr1000g_t2_avg = missing2nan(df."MR1000G-RDL_T2_AVG");
        mr1000g_sh1_std = missing2nan(df."MR1000G-RDL_SH1_STD");
        mr1000g_sh2_std = missing2nan(df."MR1000G-RDL_SH2_STD");
        mr1000g_press_avg = missing2nan(df."MR1000G-RDL_PRESS_AVG");
        mr1000g_incly_avg = missing2nan(df."MR1000G-RDL_INCLY_AVG");
        mr1000g_eps1 = missing2nan(df."MR1000G-RDL_EPS1");
        mr1000g_qc1 = missing2nan(df."MR1000G-RDL_QC1");
        mr1000g_eps2 = missing2nan(df."MR1000G-RDL_EPS2");
        mr1000g_qc2 = missing2nan(df."MR1000G-RDL_QC2");

        tmpvar = Array{Float64,1}(undef, length(ad2cp_alt));
        tmpvar .= NaN;

        if dataflag == "all"
            ad2cp_Unorth = missing2nan(df.AD2CP_Unorth_c);
            ad2cp_Ueast = missing2nan(df.AD2CP_Ueast_c);
            ad2cp_Utot = missing2nan(df.AD2CP_Utot_c);
            ad2cp_Udir = missing2nan(df.AD2CP_Udir_c);
            ad2cp_qf = missing2nan(df.AD2CP_QF_c);
        else
            ad2cp_Unorth = tmpvar;
            ad2cp_Ueast = tmpvar;
            ad2cp_Utot = tmpvar;
            ad2cp_Udir = tmpvar;
            ad2cp_qf = tmpvar;
        end

        display(yo)

        yo1d = cat(yo1d, yo, dims = 1);
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
        
        push!(pld_rt, PLD_RT(gliderSN, missionID, project, yo, t, z, lon, lat, nav_resource, ad2cp_time, ad2cp_heading, ad2cp_pitch, ad2cp_roll, ad2cp_pressure, ad2cp_alt,  ad2cp_v1_cn1, ad2cp_v2_cn1, ad2cp_v3_cn1, ad2cp_v4_cn1, ad2cp_v1_cn2, ad2cp_v2_cn2, ad2cp_v3_cn2, ad2cp_v4_cn2, ad2cp_v1_cn3, ad2cp_v2_cn3, ad2cp_v3_cn3, ad2cp_v4_cn3, ad2cp_v1_cn4, ad2cp_v2_cn4, ad2cp_v3_cn4, ad2cp_v4_cn4, ad2cp_v1_cn5, ad2cp_v2_cn5, ad2cp_v3_cn5, ad2cp_v4_cn5, ad2cp_v1_cn6, ad2cp_v2_cn6, ad2cp_v3_cn6, ad2cp_v4_cn6, flbbcd_chl_count, flbbcd_chl_scaled, flbbcd_bb_700_count, flbbcd_bb_700_scaled, flbbcd_cdom_count, flbbcd_cdom_scaled, legato_conductivity, legato_temperature, legato_pressure, legato_salinity, legato_condtemp, mr1000g_t1_avg, mr1000g_t2_avg, mr1000g_sh1_std, mr1000g_sh2_std, mr1000g_press_avg, mr1000g_incly_avg, mr1000g_eps1, mr1000g_qc1, mr1000g_eps2, mr1000g_qc2, ad2cp_Unorth, ad2cp_Ueast, ad2cp_Utot, ad2cp_Udir, ad2cp_qf));    
    end #for
    pld1d_rt = PLD_RT(gliderSN, missionID, project, yo1d, t1d, z1d, lon1d, lat1d, nav_resource1d, ad2cp_time1d, ad2cp_heading1d, ad2cp_pitch1d, ad2cp_roll1d, ad2cp_pressure1d, ad2cp_alt1d, ad2cp_v1_cn1_1d, ad2cp_v2_cn1_1d, ad2cp_v3_cn1_1d, ad2cp_v4_cn1_1d, ad2cp_v1_cn2_1d, ad2cp_v2_cn2_1d, ad2cp_v3_cn2_1d, ad2cp_v4_cn2_1d, ad2cp_v1_cn3_1d, ad2cp_v2_cn3_1d, ad2cp_v3_cn3_1d, ad2cp_v4_cn3_1d, ad2cp_v1_cn4_1d, ad2cp_v2_cn4_1d, ad2cp_v3_cn4_1d, ad2cp_v4_cn4_1d, ad2cp_v1_cn5_1d, ad2cp_v2_cn5_1d, ad2cp_v3_cn5_1d, ad2cp_v4_cn5_1d, ad2cp_v1_cn6_1d, ad2cp_v2_cn6_1d, ad2cp_v3_cn6_1d, ad2cp_v4_cn6_1d, flbbcd_chl_count_1d, flbbcd_chl_scaled_1d, flbbcd_bb_700_count_1d, flbbcd_bb_700_scaled_1d, flbbcd_cdom_count_1d, flbbcd_cdom_scaled_1d, legato_conductivity_1d, legato_temperature_1d, legato_pressure_1d, legato_salinity_1d, legato_condtemp_1d, mr1000g_t1_avg_1d, mr1000g_t2_avg_1d, mr1000g_sh1_std_1d, mr1000g_sh2_std_1d, mr1000g_press_avg_1d, mr1000g_incly_avg_1d, mr1000g_eps1_1d, mr1000g_qc1_1d, mr1000g_eps2_1d, mr1000g_qc2_1d, ad2cp_Unorth_1d, ad2cp_Ueast_1d, ad2cp_Utot_1d, ad2cp_Udir_1d, ad2cp_qf_1d);

    # combinating NAV_RT and PLD_RT data into one glider data structure
    #gliderRT = SeaExplorerRT(nav_rt, pld_rt, nav1d_rt, pld1d_rt);

    return pld_rt, pld1d_rt
end

function seaexplorer_load_mission(gliderSN, missionID, project)

    # setting src and data directory paths
    srcdir = "/Users/gong/GitHub/jlglider/";
    dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
    #dataroot = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/2022_fieldwork/";
    #dataroot = "/Users/gong/Research/sea064/";

    if (@isdefined gliderSN) != true
        gliderSN = "SEA064"
    end

    if (@isdefined missionID) != true
        mission = 37
    end

    if (@isdefined project) != true
        project = "norse-janmayen"
    end

    # define data load location
    #datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";

    # define dataset loading parameters
    #project = "maracoos"
    #deploydate = "20220311"
    #suffix = "data"

    if (gliderSN == "SEA064") & (missionID == 37)
        projectname = "norse"
        deploydate = "20221021"
        suffix = "deployment"
        datadir = dataroot * "sea064-20221021-norse-janmayen-complete/";
        dataflag = "all";
    elseif (gliderSN == "SEA064") & (missionID == 38)
        projectname = "norse"
        deploydate = "20221102"
        suffix = "deployment"
        datadir = dataroot * "sea064-20221102-norse-lofoten-complete/";
        dataflag = "all";
    elseif (gliderSN == "SEA064") & (missionID == 48)
        projectname = "norse"
        deploydate = "20231112"
        suffix = "deployment"
        datadir = dataroot * "sea064-20231112-norse/";
        dataflag = "all";
    elseif (gliderSN == "SEA064") & (missionID == 58)
        projectname = "nesma"
        deploydate = "20240709"
        suffix = "deployment"
        datadir = dataroot * "sea064-20240709-nesma";
        dataflag = "realtime";
    elseif (gliderSN == "SEA094") & (missionID == 41)
        projectname = "nesma"
        deploydate = "20240709"
        suffix = "deployment"
        datadir = dataroot * "sea094-20240709-nesma";
        dataflag = "realtime";
    end

    if (dataflag == "realtime") | (dataflag == "all")
        navdir = datadir * "glimpse/";
        scidir = datadir * "glimpse/";
    else 
        navdir = datadir * "nav/logs/";
        scidir = datadir * "science/logs/";
    end

    (SEAnav, SEAnav1d) = load_NAV(gliderSN, missionID, project, navdir, dataflag);
    (SEApld, SEApld1d) = load_PLD(gliderSN, missionID, project, scidir, dataflag); # last dataflag parameter, 0 for sub individual files, 1 for sub all, >2 for raw individual files

    return SEAnav, SEAnav1d, SEApld, SEApld1d
end

function seaexplorer_load_mission(missionYAML::String)

    mission = YAML.load_file(missionYAML);
    gliderName = mission["gliderName"];
    gliderSN = mission["gliderSN"];
    missionID = mission["missionID"];
    project = mission["project"];
    deploydate = mission["deploydate"];
    dataroot = mission["dataroot"];
    suffix = mission["suffix"];
    dataflag = mission["dataflag"];
    
    datadir = dataroot * lowercase(gliderName) * "-" * string(deploydate) * "-" * lowercase(project) * "-" * suffix * "/";

    if (dataflag == "realtime") | (dataflag == "all")
        navdir = datadir * "glimpse/";
        scidir = datadir * "glimpse/";
    else 
        navdir = datadir * "nav/logs/";
        scidir = datadir * "science/logs/";
    end

    (SEAnav, SEAnav1d) = load_NAV(gliderSN, missionID, project, navdir, dataflag);
    (SEApld, SEApld1d) = load_PLD(gliderSN, missionID, project, scidir, dataflag); # last dataflag parameter, 0 for sub individual files, 1 for sub all, >2 for raw individual files

    return SEAnav, SEAnav1d, SEApld, SEApld1d
end

function seaexplorer_process(sea064pld1d)
    lon = sea064pld1d.lon;
    lat = sea064pld1d.lat;

    gind = findall(sea064pld1d.nav_resource .== 110);
    gpst = sea064pld1d.t[gind];
    gpslon, gpslat = lon[gind], lat[gind];

    gliderSN = sea064pld1d.gliderSN;
    missionID = sea064pld1d.missionID;
    project = sea064pld1d.project;

    t = sea064pld1d.t;
    yo = sea064pld1d.yo;
    ns = sea064pld1d.nav_resource;

    badind = findall((lon .== 0.0 .&& lat .== 0.0) .|| (t .< DateTime(2020,1,1,0,0,0)));
    if isempty(badind) != true
        lon[badind] .= NaN;
        lat[badind] .= NaN;
        t[badind] .= t[badind[end]+1];
    end

    #t = cleanTime(sea064pld1d.t);
    p = cleanPress(sea064pld1d.legato_pressure);
    z = gsw_z_from_p.(p, lat, 0, 0);
    temp = cleanTemp(sea064pld1d.legato_temperature);
    salt = cleanSalt(sea064pld1d.legato_salinity);
    saltA = cleanSalt(gsw_sa_from_sp.(salt, p, lon, lat));
    ctemp = cleanTemp(gsw_ct_from_t.(saltA, temp, p));
    sigma0 = sigma0_from_t_sp(temp, salt, p, lon, lat);
    spice0 = spice0_from_t_sp(temp, salt, p, lon, lat);
    sndspd = gsw_sound_speed.(saltA, ctemp, p);

    n2, pmid = N2_from_t_sp(temp, salt, p, lon, lat);
    zmid = gsw_z_from_p.(pmid, lat[2:end], 0, 0);
    tmid = t[1:end-1] .+ Second(15);

    mr_eps1 = cleanEPS(sea064pld1d.mr1000g_eps1);
    mr_eps2 = cleanEPS(sea064pld1d.mr1000g_eps2);
    mr_qc1 = clean9999(sea064pld1d.mr1000g_qc1);
    mr_qc2 = clean9999(sea064pld1d.mr1000g_qc2);
    mr_sh1_std = clean9999(sea064pld1d.mr1000g_sh1_std);
    mr_sh2_std = clean9999(sea064pld1d.mr1000g_sh2_std);
    mr_t1_avg = clean9999(sea064pld1d.mr1000g_t1_avg);
    mr_t2_avg = clean9999(sea064pld1d.mr1000g_t2_avg);
    ad2cp_Ueast = cleanAD2CP(sea064pld1d.ad2cp_Ueast);
    ad2cp_Unorth = cleanAD2CP(sea064pld1d.ad2cp_Unorth);
    ad2cp_Utot = cleanAD2CP(sea064pld1d.ad2cp_Utot);
    ad2cp_Udir = cleanAD2CP(sea064pld1d.ad2cp_Udir);
    ad2cp_qf = cleanAD2CP(sea064pld1d.ad2cp_qf);
    chla = cleanFLBBCDchl(sea064pld1d.flbbcd_chl_scaled);
    bb700 = cleanFLBBCDbb700(sea064pld1d.flbbcd_bb_700_scaled);
    cdom = cleanFLBBCDcdom(sea064pld1d.flbbcd_cdom_scaled);

    SEAdata = SeaExplorerData(gliderSN, missionID, project, yo, ns, t, lon, lat, gpst, gpslon, gpslat, p, z, temp, salt, saltA, ctemp, sigma0, spice0, sndspd, mr_eps1, mr_eps2, mr_qc1, mr_qc2, mr_sh1_std, mr_sh2_std, mr_t1_avg, mr_t2_avg, ad2cp_Ueast, ad2cp_Unorth, ad2cp_Utot, ad2cp_Udir, ad2cp_qf, chla, bb700, cdom, n2, pmid, zmid, tmid);
end


#=
function seaexplorer_process(sea064pld1d)
    lon = sea064pld1d.lon;
    lat = sea064pld1d.lat;
    t = sea064pld1d.t;
    yo = sea064pld1d.yo;
    ns = sea064pld1d.nav_resource;

    badind = findall((lon .== 0.0 .&& lat .== 0.0) .|| (t .< DateTime(2020,1,1,0,0,0)));
    if isempty(badind) != true
        lon[badind] .= NaN;
        lat[badind] .= NaN;
        t[badind] .= t[badind[end]+1];
    end

    #t = cleanTime(sea064pld1d.t);
    p = cleanPress(sea064pld1d.legato_pressure);
    z = gsw_z_from_p.(p, lat, 0, 0);
    temp = cleanTemp(sea064pld1d.legato_temperature);
    salt = cleanSalt(sea064pld1d.legato_salinity);
    saltA = cleanSalt(gsw_sa_from_sp.(salt, p, lon, lat));
    ctemp = cleanTemp(gsw_ct_from_t.(saltA, temp, p));
    sigma0 = sigma0_from_t_sp(temp, salt, p, lon, lat);
    spice0 = spice0_from_t_sp(temp, salt, p, lon, lat);
    sndspd = gsw_sound_speed.(saltA, ctemp, p);

    n2, pmid = N2_from_t_sp(temp, salt, p, lon, lat);
    zmid = gsw_z_from_p.(pmid, lat[2:end], 0, 0);
    tmid = t[1:end-1] .+ Second(15);

    mr_eps1 = cleanEPS(sea064pld1d.mr1000g_eps1);
    mr_eps2 = cleanEPS(sea064pld1d.mr1000g_eps2);
    mr_qc1 = clean9999(sea064pld1d.mr1000g_qc1);
    mr_qc2 = clean9999(sea064pld1d.mr1000g_qc2);
    mr_sh1_std = clean9999(sea064pld1d.mr1000g_sh1_std);
    mr_sh2_std = clean9999(sea064pld1d.mr1000g_sh2_std);
    mr_t1_avg = clean9999(sea064pld1d.mr1000g_t1_avg);
    mr_t2_avg = clean9999(sea064pld1d.mr1000g_t2_avg);
    ad2cp_Ueast = cleanAD2CP(sea064pld1d.ad2cp_Ueast);
    ad2cp_Unorth = cleanAD2CP(sea064pld1d.ad2cp_Unorth);
    ad2cp_Utot = cleanAD2CP(sea064pld1d.ad2cp_Utot);
    ad2cp_Udir = cleanAD2CP(sea064pld1d.ad2cp_Udir);
    ad2cp_qf = cleanAD2CP(sea064pld1d.ad2cp_qf);
    chla = cleanFLBBCDchl(sea064pld1d.flbbcd_chl_scaled);
    bb700 = cleanFLBBCDbb700(sea064pld1d.flbbcd_bb_700_scaled);
    cdom = cleanFLBBCDcdom(sea064pld1d.flbbcd_cdom_scaled);

    SEAdata = SeaExplorerData(yo, t, ns, lon, lat, p, z, temp, salt, saltA, ctemp, sigma0, spice0, sndspd, mr_eps1, mr_eps2, mr_qc1, mr_qc2, mr_sh1_std, mr_sh2_std, mr_t1_avg, mr_t2_avg, ad2cp_Ueast, ad2cp_Unorth, ad2cp_Utot, ad2cp_Udir, ad2cp_qf, chla, bb700, cdom, n2, pmid, zmid, tmid);
end
=#

function seaexplorer_MR_laur_load(pld1, pld2, pz, dz)
    # load roughly processed MR1000G data by Laur Ferris for plotting
    # gong@vims.edu 2023-04-28

    #if (@isdefined jmpld1d) != true
    #    include("run_seaexplorer.jl");
    #end
    
    #pld1 = jmpld1d;
    #pld2 = lbepld1d;
    
    pld1t = datetime2unix.(pld1.t);
    pld2t = datetime2unix.(pld2.t);
    
    # load performatted epsilon values produced by Laur Ferris
    datadir = "/Users/gong/GitHub/jlglider/seaexplorer/data/"
    mr1 = matopen(datadir * "jm_P.mat");
    mr2 = matopen(datadir * "lb_P.mat");
    
    mr1_t = read(mr1, "unixt");
    mr1_t[:,106] .= datetime2unix(DateTime(2022,10,30,1,0,0));
    mr1_p = read(mr1, "pres");
    mr1_eps1 = read(mr1, "eps1");
    mr1_eps2 = read(mr1, "eps2");
    
    mr2_t = read(mr2, "unixt");
    mr2_p = read(mr2, "pres");
    mr2_eps1 = read(mr2, "eps1");
    mr2_eps2 = read(mr2, "eps2");
    
    # calculate z from p assuming latitude of 70
    mr1_z = gsw_z_from_p.(mr1_p, 70, 0, 0);
    mr2_z = gsw_z_from_p.(mr2_p, 70, 0, 0);
    
    # create interpolation function for lon and lat based on time using glider CTD data
    lon1f = linear_interpolation(pld1t, pld1.lon);
    lat1f = linear_interpolation(pld1t, pld1.lat);
    lon2f = linear_interpolation(pld2t, pld2.lon, extrapolation_bc=Line());
    lat2f = linear_interpolation(pld2t, pld2.lat, extrapolation_bc=Line()); # extrapolation_bc=Line()

    # calculate lon and lat for the epsilon casts based on interpolation functions
    lon1 = zeros(Float64, size(mr1_t,1), size(mr1_t,2));
    lat1 = deepcopy(lon1);
    lon2 = zeros(Float64, size(mr2_t,1), size(mr2_t,2));
    lat2 = deepcopy(lon2);
    
    eps1 = zeros(Float64, size(mr1_t,2));
    eps2 = zeros(Float64, size(mr2_t,2));
    
    # calculate depth-averaged epsilon values based on depth range of zlo and zhi
    zlo, zhi = pz-dz, pz+dz;

    if abs(zlo) < 0.1
        zlo = -0.1;
    end

    for i = 1:size(mr1_p, 2)
        lon1[:,i] = lon1f(mr1_t[:,i]);
        lat1[:,i] = lat1f(mr1_t[:,i]);
        bind = findall(abs.(mr1_z[:,i]) .< 0.01);
        mr1_eps1[bind,i] .= NaN;
        mr1_eps2[bind,i] .= NaN;
    
        zind = findall(zlo .<= mr1_z[:,i] .<= zhi);
        #eps1[i] = 10 .^ NaNMath.mean(log10.(mr1_eps1[zind,i]));
        eps1[i] = 10 .^ NaNMath.mean([NaNMath.mean(log10.(mr1_eps1[zind,i])); NaNMath.mean(log10.(mr1_eps2[zind,i]))]);
    end
    
    for i = 1:size(mr2_p, 2)
        lon2[:,i] = lon2f(mr2_t[:,i]);
        lat2[:,i] = lat2f(mr2_t[:,i]);
        bind = findall(abs.(mr2_z[:,i]) .< 0.01);
        mr2_eps1[bind,i] .= NaN;
        mr2_eps2[bind,i] .= NaN;
   
        zind = findall(zlo .<= mr2_z[:,i] .<= zhi);
        #eps2[i] = 10 .^ NaNMath.mean(log10.(mr2_eps1[zind,i]));
        eps2[i] = 10 .^ NaNMath.mean([NaNMath.mean(log10.(mr2_eps1[zind,i])); NaNMath.mean(log10.(mr2_eps2[zind,i]))]);    
    end
    
    return lon1, lat1, lon2, lat2, eps1, eps2
end    

end