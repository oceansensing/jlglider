using Glob, DataFrames, CSV, Dates

mutable struct NAV
    time::Array{DateTime};
    depth::Array{AbstractFloat};
    lon::Array{AbstractFloat};
    lat::Array{AbstractFloat};
end

mutable struct LEGATO
    time::Array{DateTime};
    pres::Array{AbstractFloat};
    temp::Array{AbstractFloat};
    cond::Array{AbstractFloat};
    condtemp::Array{AbstractFloat};
    salt::Array{AbstractFloat};
end

mutable struct FLBBCD
    time::Array{DateTime};
    chla::Array{AbstractFloat};
    cdom::Array{AbstractFloat};
    bb700::Array{AbstractFloat};
end

mutable struct SeaExplorer
    nav::NAV;
    ctd::LEGATO;
    flbbcd::FLBBCD;
end

srcdir = "/Users/gong/Research/jlglider/";
dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";

project = "maracoos"
deploydate = "20220311"
suffix = "data"
glidername = "sea064"
mission = "24"

missionroot = glidername * "." * mission;
gliroot = missionroot * "." * "gli.sub.";
pldroot_rt = missionroot * "." * "pld1.sub.";
pldroot_raw = missionroot * "." * "pld1.raw.";
ad2cproot_raw = missionroot * "." * "ad2cp.raw.";
legatoroot_raw = missionroot * "." * "legato.raw.";

datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
navdir = datadir * "nav/";
scidir = datadir * "science/";

glilist = Glob.glob(gliroot * "*", navdir);
pldlist_rt = Glob.glob(pldroot_rt * "*", scidir);
pldlist_raw = Glob.glob(pldroot_raw * "*", scidir);
ad2cplist_raw = Glob.glob(ad2cproot_raw * "*", scidir);
legatolist_raw = Glob.glob(legatoroot_raw * "*", scidir);

timeformat = "dd/mm/yyyy HH:MM:SS.sss"
global time1d = [];
global lon1d = [];
global lat1d = [];

#i = 1
#print(scidir * pldroot_raw * string(i) * "\n")
#df = CSV.read(scidir * pldroot_raw * string(i), header=1, delim=";", DataFrame);
#global time1d = cat(time1d, DateTime.(df.PLD_REALTIMECLOCK, timeformat), dims=1);

for i = 1:length(pldlist_raw)
#for i = 1:1
        display(i)
    print(scidir * pldroot_raw * string(i) * "\n")
    df = CSV.read(scidir * pldroot_raw * string(i), header=1, delim=";", DataFrame);

    navlon = df.NAV_LONGITUDE;
    navlat = df.NAV_LATITUDE;
    lon = Array{Float64,1};
    lat = Array{Float64,1};
    lon = trunc.(navlon ./ 100) + (navlon .- trunc.(navlon ./ 100)*100) / 60;
    lat = trunc.(navlat ./ 100) + (navlat .- trunc.(navlat ./ 100)*100) / 60;
    
    global time1d = cat(time1d, DateTime.(df.PLD_REALTIMECLOCK, timeformat), dims=1);
    global lon1d = cat(lon1d, lon, dims=1);
    global lat1d = cat(lat1d, lat, dims=1);
end
