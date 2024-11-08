# seaexplorerYAMLtest.jl
# DG 2024-11-07

include("seaexplorerType.jl");
include("gsw_c2po.jl");
include("seaexplorerFunc.jl");

using Glob, YAML
using Glob, DataFrames, CSV, Dates, Missings, NaNMath, Interpolations, YAML
using NCDatasets
using GibbsSeaWater, MAT
import .seaexplorerType: NAV_RT, PLD_RT, SeaExplorerData, SeaExplorerCTD, LEGATO
import .gsw_c2po: sigma0_from_t_sp, spice0_from_t_sp, N2_from_t_sp
import .seaexplorerFunc: load_NAV, load_PLD, load_LEGATO, seaexplorer_load_mission, seaexplorer_process, missing2nan
import .seaexplorerFunc: cleanTime, cleanPress, cleanTemp, clean9999, cleanAD2CP, cleanFLBBCDchl, cleanFLBBCDbb700, cleanFLBBCDcdom, cleanSalt, cleanEPS

i = 2

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdirpath = "/Users/gong/GitHub/jlglider/seaexplorer/mission_yaml_PASSENGERS/";

if missionYAMLdirpath[end-4:end] != ".yaml"
    if missionYAMLdirpath[end] != '/'
     missionYAMLdirpath *= "/";
    end
    missionYAMLpath = Glob.glob("*.yaml", missionYAMLdirpath)
else
    missionYAMLpath = [missionYAMLdirpath];
end    

#gliderCTDarray = SeaExplorerData[];

display(missionYAMLpath[i])
#SEAnav, SEAnav1d, SEApld, SEApld1d = seaexplorer_load_mission(missionYAMLpath[i]);
## BEGIN seaexplorer_load_mission
mission = YAML.load_file(missionYAMLpath[i]);
gliderName = mission["gliderName"];
gliderSN = mission["gliderSN"];
missionID = mission["missionID"];
project = mission["project"];
deploydate = mission["deploydate"];
dataroot = mission["dataroot"];
suffix = mission["suffix"];
dataflag = mission["dataflag"];

datadir = dataroot * lowercase(gliderName) * "-" * string(deploydate) * "-" * lowercase(project) * "-" * suffix * "/";

if (dataflag == "realtime") || (dataflag == "all")
    navdir = datadir * "glimpse/";
    scidir = datadir * "glimpse/";
else 
    navdir = datadir * "delayed/nav/logs/";
    if (gliderSN == 94) && (missionID == 41)
        scidir = datadir * "delayed/pld064/home/user/logs/";
    else
        scidir = datadir * "delayed/pld1/logs/";
    end
end

#(SEAnav, SEAnav1d) = load_NAV(gliderSN, missionID, project, navdir, dataflag);
## BEGIN load_NAV
SEAnav = NAV_RT[];

#missionroot = uppercase(glidername) * "." * mission;
#gliroot_rt = missionroot * "." * "gli.sub.";
#glilist_rt = Glob.glob(gliroot_rt * "*", navdir);

global glilist = Glob.glob( "*" * string(gliderSN; pad=3) * "." * string(missionID) * ".gli.sub.*", navdir);

# separating '.all' from '.###' files
#glilist_suffix = [];
global yos = [];
global yolist = [];
allindx = 1;
for i = 1:length(glilist)
    suffix2 = glilist[i][end-1:end]
    fnlen = length(glilist[i]) - length(navdir);
    if suffix2 == "gz"
        if fnlen == 22  # sea064.37.gli.sub.1.gz
            global yo = parse(Int,glilist[i][end-3:end-3]);
        elseif fnlen == 23 # sea064.37.gli.sub.10.gz
            global yo = parse(Int,glilist[i][end-4:end-3]);
        elseif fnlen == 24 # sea064.37.gli.sub.100.gz
            global yo = parse(Int,glilist[i][end-5:end-3]);
        elseif fnlen == 25 # sea064.37.gli.sub.1000.gz
            global yo = parse(Int,glilist[i][end-6:end-3]);
        end
        global yos = push!(yos, yo);
        global yolist = push!(yolist, glilist[i]);
    elseif (suffix2 != "gz") & (suffix2 != "ll") & (suffix2 != "sv")
        if fnlen == 19  # sea064.37.gli.sub.1
            global yo = parse(Int,glilist[i][end:end]);
        elseif fnlen == 20 # sea064.37.gli.sub.10
            global yo = parse(Int,glilist[i][end-1:end]);
        elseif fnlen == 21 # sea064.37.gli.sub.100
            global yo = parse(Int,glilist[i][end-2:end]);
        elseif fnlen == 22 # sea064.37.gli.sub.1000
            global yo = parse(Int,glilist[i][end-3:end]);
        end
        global yos = push!(yos, yo);
        global yolist = push!(yolist, glilist[i]);
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
global yo1d = [];
global t1d = [];
global z1d = [];
global lon1d = [];
global lat1d = [];
global NavState1d = [];
global SecurityLevel1d = [];
global Heading1d = [];
global Declination1d = [];
global Pitch1d = [];
global Roll1d = [];
global DeadReckoning1d = [];
global DesiredH1d = [];
global BallastCmd1d = [];
global BallastPos1d = [];
global LinCmd1d = [];
global LinPos1d = [];
global AngCmd1d = [];
global AngPos1d = [];
global Voltage1d = [];
global Altitude1d = [];    

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

    global yo1d = cat(yo1d, yo, dims=1);
    global t1d = cat(t1d, t, dims=1);
    global z1d = cat(z1d, z, dims=1);
    global lon1d = cat(lon1d, lon, dims=1);
    global lat1d = cat(lat1d, lat, dims=1);
    global NavState1d = cat(NavState1d, NavState, dims=1);
    global SecurityLevel1d = cat(SecurityLevel1d, SecurityLevel, dims=1);
    global Heading1d = cat(Heading1d, Heading, dims=1);
    global Declination1d = cat(Declination1d, Declination, dims=1);
    global Pitch1d = cat(Pitch1d, Pitch, dims=1);
    global Roll1d = cat(Roll1d, Roll, dims=1);
    global DeadReckoning1d = cat(DeadReckoning1d, DeadReckoning, dims=1);
    global DesiredH1d = cat(DesiredH1d, DesiredH, dims=1);
    global BallastCmd1d = cat(BallastCmd1d, BallastCmd, dims=1);
    global BallastPos1d = cat(BallastPos1d, BallastPos, dims=1);
    global LinCmd1d = cat(LinCmd1d, LinCmd, dims=1);
    global LinPos1d = cat(LinPos1d, LinPos, dims=1);
    global AngCmd1d = cat(AngCmd1d, AngCmd, dims=1);
    global AngPos1d = cat(AngPos1d, AngPos, dims=1);
    global Voltage1d = cat(Voltage1d, Voltage, dims=1);
    global Altitude1d = cat(Altitude1d, Altitude, dims=1);

    push!(SEAnav, NAV_RT(gliderSN, missionID, project, yo, t, z, lon, lat, NavState, SecurityLevel, Heading, Declination, Pitch, Roll, DeadReckoning, DesiredH, BallastCmd, BallastPos, LinCmd, LinPos, AngCmd, AngPos, Voltage, Altitude));
end #for
SEAnav1d = NAV_RT(gliderSN, missionID, project, yo1d, t1d, z1d, lon1d, lat1d, NavState1d, SecurityLevel1d, Heading1d, Declination1d, Pitch1d, Roll1d, DeadReckoning1d, DesiredH1d, BallastCmd1d, BallastPos1d, LinCmd1d, LinPos1d, AngCmd1d, AngPos1d, Voltage1d, Altitude1d);
## END load_NAV

#(SEApld, SEApld1d) = load_LEGATO(gliderSN, missionID, project, scidir, dataflag); # last dataflag parameter, 0 for sub individual files, 1 for sub all, >2 for raw individual files
## END seaexplorer_load_mission

## BEGIN load_LEGATO
legato_flag = 1;
SEApld = LEGATO[];
if (dataflag == "sub") | (dataflag == "realtime") | (dataflag == "all")
    datatype = "sub"
elseif (dataflag == "raw") | (dataflag == "delayed")
    datatype = "raw"
end
pldlist = Glob.glob("*" * string(gliderSN; pad=3) * "." * string(missionID) * ".pld1." * datatype * ".*", scidir);

global yos = [];
global yolist = [];
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
        global yos = push!(yos, yo);
        global yolist = push!(yolist, pldlist[i]);
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
        global yos = push!(yos, yo);
        global yolist = push!(yolist, pldlist[i]);
    elseif suffix2 == "ll"
        if fnlen == 21+1 # SEA064.37.pld1.sub.all
            allindx = i;
        end
        global yos = [-1];
    elseif suffix2 == "sv"
        if fnlen == 25 # SEA064.37.gli.sub.all.csv
            allindx = i;
        end
        global yos = [-1];
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

timeformat = "dd/mm/yyyy HH:MM:SS.sss"

global yo1d = [];
global t1d = [];
global z1d = [];
global lon1d = [];
global lat1d = [];
global nav_resource1d = [];
global legato_conductivity_1d = [];
global legato_temperature_1d = [];
global legato_pressure_1d = [];
global legato_salinity_1d = [];
global legato_condtemp_1d = [];

for i = 1:length(pldlist)
#i = 56
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
    nav_resource = missing2nan(df.NAV_RESOURCE);

    global yo1d = cat(yo1d, yo, dims = 1);
    global t1d = cat(t1d, t, dims = 1);
    global z1d = cat(z1d, z, dims = 1);
    global lon1d = cat(lon1d, lon, dims = 1);
    global lat1d = cat(lat1d, lat, dims = 1);
    global nav_resource1d = cat(nav_resource1d, nav_resource, dims = 1);

    if legato_flag == 1
        legato_conductivity = missing2nan(df.LEGATO_CONDUCTIVITY);
        legato_temperature = missing2nan(df.LEGATO_TEMPERATURE);
        legato_pressure = missing2nan(df.LEGATO_PRESSURE);
        legato_salinity = missing2nan(df.LEGATO_SALINITY);
        legato_condtemp = missing2nan(df.LEGATO_CONDTEMP);

        global legato_conductivity_1d = cat(legato_conductivity_1d, legato_conductivity, dims = 1);
        global legato_temperature_1d = cat(legato_temperature_1d, legato_temperature, dims = 1);
        global legato_pressure_1d = cat(legato_pressure_1d, legato_pressure, dims = 1);
        global legato_salinity_1d = cat(legato_salinity_1d, legato_salinity, dims = 1);
        global legato_condtemp_1d = cat(legato_condtemp_1d, legato_condtemp, dims = 1);    
    end   

    display(yo)
    
    push!(SEApld, LEGATO(gliderSN, missionID, project, yo, t, z, lon, lat, nav_resource, legato_conductivity, legato_temperature, legato_pressure, legato_salinity, legato_condtemp));    
end #for
SEApld1d = LEGATO(gliderSN, missionID, project, yo1d, t1d, z1d, lon1d, lat1d, nav_resource1d, legato_conductivity_1d, legato_temperature_1d, legato_pressure_1d, legato_salinity_1d, legato_condtemp_1d);

## END load_LEGATO

## BEGIN seaexplorer_process
lon = SEApld1d.lon;
lat = SEApld1d.lat;

gind = findall(SEApld1d.nav_resource .== 110);
gpst = SEApld1d.t[gind];
gpslon, gpslat = lon[gind], lat[gind];

gliderSN = SEApld1d.gliderSN;
missionID = SEApld1d.missionID;
project = SEApld1d.project;

t = SEApld1d.t;
yo = SEApld1d.yo;
ns = SEApld1d.nav_resource;

badind = findall((lon .== 0.0 .&& lat .== 0.0) .|| (t .< DateTime(2020,1,1,0,0,0)));
if isempty(badind) != true
    lon[badind] .= NaN;
    lat[badind] .= NaN;
    t[badind] .= t[badind[end]+1];
end

#t = cleanTime(SEApld1d.t);
p = cleanPress(SEApld1d.legato_pressure);
z = gsw_z_from_p.(p, lat, 0, 0);
temp = cleanTemp(SEApld1d.legato_temperature);
salt = cleanSalt(SEApld1d.legato_salinity);
saltA = cleanSalt(gsw_sa_from_sp.(salt, p, lon, lat));
ctemp = cleanTemp(gsw_ct_from_t.(saltA, temp, p));
sigma0 = sigma0_from_t_sp(temp, salt, p, lon, lat);
spice0 = spice0_from_t_sp(temp, salt, p, lon, lat);
sndspd = gsw_sound_speed.(saltA, ctemp, p);

n2, pmid = N2_from_t_sp(temp, salt, p, lon, lat);
zmid = gsw_z_from_p.(pmid, lat[2:end], 0, 0);
tmid = t[1:end-1] .+ Second(15);

SEAdata = SeaExplorerCTD(gliderSN, missionID, project, yo, ns, t, lon, lat, gpst, gpslon, gpslat, p, z, temp, salt, saltA, ctemp, sigma0, spice0, sndspd, n2, pmid, zmid, tmid);
## END load_LEGATO