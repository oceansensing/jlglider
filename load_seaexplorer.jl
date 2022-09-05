# This script load SeaExplorer data file
# 2022-07-16: Donglai Gong

using Glob, DataFrames, CSV, Dates, Missings

# creating glider data types
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
    alt::Array{AbstractFloat};
    heading::Array{AbstractFloat};
    pitch::Array{AbstractFloat};
    pressure::Array{AbstractFloat};
    roll::Array{AbstractFloat};
    v1_cn2::Array{AbstractFloat};
    v1_cn3::Array{AbstractFloat};
end

mutable struct SeaExplorer
    nav::Array{NAV};
    ctd::Array{LEGATO};
    flbbcd::Array{FLBBCD};
    nav1d::NAV;
    ctd1d::LEGATO;
    flbbcd1d::FLBBCD;
end

# define a function that converts missings in an array to NaN
function missing2nan(varin)
    varout = Float64.(collect(Missings.replace(varin, NaN)));
end

function load_SEdata(glidername::String, mission::String, navdir::String, scidir::String)

    # initilized data strucctures
    nav = NAV[];
    ctd = LEGATO[];
    flbbcd = FLBBCD[];
    ad2cp = AD2CP[];
    #glider = SeaExplorer[];

    # set directory paths
    missionroot = glidername * "." * mission;
    gliroot = missionroot * "." * "gli.sub.";
    pldroot_rt = missionroot * "." * "pld1.sub.";
    pldroot_raw = missionroot * "." * "pld1.raw.";
    ad2cproot_raw = missionroot * "." * "ad2cp.raw.";
    legatoroot_raw = missionroot * "." * "legato.raw.";

    # load data file lists
    glilist = Glob.glob(gliroot * "*", navdir);
    pldlist_rt = Glob.glob(pldroot_rt * "*", scidir);
    pldlist_raw = Glob.glob(pldroot_raw * "*", scidir);
    ad2cplist_raw = Glob.glob(ad2cproot_raw * "*", scidir);
    legatolist_raw = Glob.glob(legatoroot_raw * "*", scidir);

    # initiate variables 
    t1d = [];
    lon1d = [];
    lat1d = [];
    z1d = [];

    temp1d = [];
    condtemp1d = [];
    cond1d = [];
    salt1d = [];
    p1d = [];

    chla1d = [];
    cdom1d = [];
    bb1d = [];

    timeformat = "dd/mm/yyyy HH:MM:SS.sss"

    # loop through the list of data files
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

        # create profile data structure (organized by yo's) for NAV, CTD, and FLBBCD data
        push!(nav, NAV(t, z, lon, lat));
        push!(ctd, LEGATO(t, p, temp, cond, condtemp, salt));
        push!(flbbcd, FLBBCD(t, chla, cdom, bb700));
    end

    # storing 1d data into NAV, CTD, and FLBBCD data structures
    nav1d = NAV(t1d, z1d, lon1d, lat1d);
    ctd1d = LEGATO(t1d, p1d, temp1d, cond1d, condtemp1d, salt1d);
    flbbcd1d = FLBBCD(t1d, chla1d, cdom1d, bb1d);

    # combinating NAV, CTD, and FLBBCD data into one glider data structure
    glider = SeaExplorer(nav, ctd, flbbcd, nav1d, ctd1d, flbbcd1d);

    return glider
end

# setting src and data directory paths
srcdir = "/Users/gong/GitHub/jlglider/";
dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";

# define dataset loading parameters
project = "maracoos"
deploydate = "20220311"
suffix = "data"

glidername = "sea064"
mission = "24"

# define data load location
datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
navdir = datadir * "nav/";
scidir = datadir * "science/";

glider = load_SEdata(glidername, mission, navdir, scidir);
