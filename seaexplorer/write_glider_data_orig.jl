# Usage Example:
# write_glider_data(glider = "SEA064", mission = 48, csvdir = "./")

workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

using Glider
import Glider.seaexplorerType: LEGATO, SeaExplorerCTD

include("seaexplorerFunc.jl")
import .seaexplorerFunc: seaexplorer_load_mission, seaexplorer_process
using Dates, DataFrames, CSV, Interpolations

function write_glider_data(; glider = "SEA064", mission = "48", csvdir = "./")
# DG: this function writes out the seaexplorer glider mission data for use with gliderad2cp library

    nav, nav1d, pld, pld1d = seaexplorer_load_mission(mission);
    seadata = seaexplorer_process(pld1d);

    nav_unixt = Dates.datetime2unix.(nav1d.t); 
    tind = sortperm(nav_unixt);
    nav_unixt = nav_unixt[tind];
    nav_dec = nav1d.Declination[tind];
    DeclinationFunc = linear_interpolation(nav_unixt, nav_dec, extrapolation_bc=Line());
    declination = DeclinationFunc(Dates.datetime2unix.(seadata.t));

    seadf = DataFrame((time = seadata.t, temperature = seadata.temp, salinity = seadata.salt, latitude = seadata.lat, pressure = seadata.p, longitude = seadata.lon, profile_number = seadata.yo, declination = declination));

    CSV.write(csvdir * glider * "_M" * string(mission) * ".csv", seadf);

end

function write_glider_data(missionYAMLdirpath::String)
    # DG: this function writes out the seaexplorer glider mission data for use with gliderad2cp library

    SEAnav, SEAnav1d, SEApld, SEApld1d = seaexplorer_load_mission(missionYAMLdirpath);
    SEApld1d = seaexplorer_process(SEApld1d);

    glidertype = SEApld1d.glidertype;
    gliderSN = SEApld1d.gliderSN;
    glidername = SEApld1d.glidername;
    missionID = SEApld1d.missionID;
    project = SEApld1d.project;

    #nav_unixt = Dates.datetime2unix.(SEAnav1d.t); 
    nav_unixt = SEAnav1d.t;
    tind = sortperm(nav_unixt);
    nav_unixt = nav_unixt[tind];
    nav_dec = SEAnav1d.Declination[tind];
    DeclinationFunc = linear_interpolation(nav_unixt, nav_dec, extrapolation_bc=Line());
    declination = DeclinationFunc(SEApld1d.t);

    seadf = DataFrame((time = Dates.unix2datetime.(SEApld1d.t), temperature = SEApld1d.temp, salinity = SEApld1d.salt, latitude = SEApld1d.lat, pressure = SEApld1d.p, longitude = SEApld1d.lon, profile_number = SEApld1d.yo, declination = declination));
    csvdir = "./"
    CSV.write(csvdir * glidername * "_M" * string(missionID) * ".csv", seadf);
end

missionYAMLdirpath = "/Users/gong/GitHub/jlglider/seaexplorer/seaexplorer_yaml_PASSENGERS/sea064-20240720-passengers.yaml";
write_glider_data(missionYAMLdirpath)