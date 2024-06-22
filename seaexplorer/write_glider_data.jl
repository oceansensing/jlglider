# Usage Example:
# write_glider_data(glider = "SEA064", mission = 48, csvdir = "./")

workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

include("seaexplorerFunc.jl")
import .seaexplorerFunc: load_NAV, load_PLD, seaexplorer_load_mission, seaexplorer_process
using Dates, DataFrames, CSV, Plots, Interpolations

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
