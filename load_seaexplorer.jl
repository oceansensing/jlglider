# This script load SeaExplorer data file
# 2022-07-16: Donglai Gong

using Glob, DataFrames, CSV, Dates, Missings

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

mutable struct SeaExplorer
    nav::Array{NAV};
    ctd::Array{LEGATO};
    flbbcd::Array{FLBBCD};
end

function missing2nan(varin)
    varout = Float64.(collect(Missings.replace(varin, NaN)));
end

nav = NAV[];
ctd = LEGATO[];
flbbcd = FLBBCD[];
#glider = SeaExplorer[];

srcdir = "/Users/gong/Research/jlglider/";
dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";

# define dataset loading parameters
project = "maracoos"
deploydate = "20220311"
suffix = "data"
glidername = "sea064"
mission = "24"

# set directory paths
missionroot = glidername * "." * mission;
gliroot = missionroot * "." * "gli.sub.";
pldroot_rt = missionroot * "." * "pld1.sub.";
pldroot_raw = missionroot * "." * "pld1.raw.";
ad2cproot_raw = missionroot * "." * "ad2cp.raw.";
legatoroot_raw = missionroot * "." * "legato.raw.";

datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
navdir = datadir * "nav/";
scidir = datadir * "science/";

# load data file lists
glilist = Glob.glob(gliroot * "*", navdir);
pldlist_rt = Glob.glob(pldroot_rt * "*", scidir);
pldlist_raw = Glob.glob(pldroot_raw * "*", scidir);
ad2cplist_raw = Glob.glob(ad2cproot_raw * "*", scidir);
legatolist_raw = Glob.glob(legatoroot_raw * "*", scidir);

timeformat = "dd/mm/yyyy HH:MM:SS.sss"
global time1d = [];
global lon1d = [];
global lat1d = [];
global z1d = [];

global temp1d = [];
global condtemp1d = [];
global cond1d = [];
global salt1d = [];
global p1d = [];

global chla1d = [];
global cdom1d = [];
global bb1d = [];

#i = 1
#print(scidir * pldroot_raw * string(i) * "\n")
#df = CSV.read(scidir * pldroot_raw * string(i), header=1, delim=";", DataFrame);
#global time1d = cat(time1d, DateTime.(df.PLD_REALTIMECLOCK, timeformat), dims=1);

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
    global time1d = cat(time1d, t, dims=1);
    global lon1d = cat(lon1d, lon, dims=1);
    global lat1d = cat(lat1d, lat, dims=1);
    global z1d = cat(z1d, z, dims=1);

    # concatinate LEGATO data into 1D arrays
    p = missing2nan(df.LEGATO_PRESSURE);
    temp = missing2nan(df.LEGATO_TEMPERATURE);
    condtemp = missing2nan(df.LEGATO_CONDTEMP);
    cond = missing2nan(df.LEGATO_CONDUCTIVITY);
    salt = missing2nan(df.LEGATO_SALINITY);

    global p1d = cat(p1d, p, dims=1);
    global temp1d = cat(temp1d, temp, dims=1);
    global condtemp1d = cat(condtemp1d, condtemp, dims=1);
    global cond1d = cat(cond1d, cond, dims=1);
    global salt1d = cat(salt1d, salt, dims=1);

    chla = missing2nan(df.FLBBCD_CHL_SCALED);
    cdom = missing2nan(df.FLBBCD_CDOM_SCALED);
    bb700 = missing2nan(df.FLBBCD_BB_700_SCALED);

    global chla1d = cat(chla1d, chla, dims=1);
    global cdom1d = cat(cdom1d, cdom, dims=1);
    global bb1d = cat(bb1d, bb700, dims=1);

    push!(nav, NAV(t, z, lon, lat));
    push!(ctd, LEGATO(t, p, temp, cond, condtemp, salt));
    push!(flbbcd, FLBBCD(t, chla, cdom, bb700));
end

glider = SeaExplorer(nav, ctd, flbbcd);

#bzind = findall(ismissing.(z1d) .== true);
#z1d[bzind] .= NaN;
#z1d = Float64.(collect(z1d));

#lon1d = Float64.(lon1d);
#lat1d = Float64.(lat1d);
#z1d = Float64.(collect(Missings.replace(z1d, NaN)));

#p1d = Float64.(collect(Missings.replace(p1d, NaN)));
#temp1d = Float64.(collect(Missings.replace(temp1d, NaN)));
#condtemp1d = Float64.(collect(Missings.replace(condtemp1d, NaN)));
#cond1d = Float64.(collect(Missings.replace(cond1d, NaN)));
#salt1d = Float64.(collect(Missings.replace(salt1d, NaN)));

