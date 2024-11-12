# Usage Example:
# write_glider_data(glider = "SEA064", mission = 48, csvdir = "./")

workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

using Glider
import Glider.seaexplorerType: LEGATO, SeaExplorerCTD

include("seaexplorerFunc.jl")
import .seaexplorerFunc: seaexplorer_load_mission, seaexplorer_process, seaexplorerCSVwrite
using Dates, DataFrames, CSV, Interpolations

# DG: this function writes out the seaexplorer glider mission data for use with gliderad2cp library
#function write_glider_data(missionYAMLdirpath::String)
#end

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLpath = [];
push!(missionYAMLpath, "/Users/gong/GitHub/jlglider/seaexplorer/seaexplorer_yaml_PASSENGERS/sea064-20240709-passengers.yaml");
push!(missionYAMLpath, "/Users/gong/GitHub/jlglider/seaexplorer/seaexplorer_yaml_PASSENGERS/sea064-20240720-passengers.yaml");
push!(missionYAMLpath, "/Users/gong/GitHub/jlglider/seaexplorer/seaexplorer_yaml_NORSE/sea064-20221021-norse.yaml");
push!(missionYAMLpath, "/Users/gong/GitHub/jlglider/seaexplorer/seaexplorer_yaml_NORSE/sea064-20221102-norse.yaml");
push!(missionYAMLpath, "/Users/gong/GitHub/jlglider/seaexplorer/seaexplorer_yaml_NORSE/sea064-20231112-norse.yaml");

seaexplorerCSVwrite(missionYAMLpath, gliderdatadir * "CSV/");

#=
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

PitchFunc = linear_interpolation(nav_unixt[tind], SEAnav1d.Pitch[tind], extrapolation_bc=Line()); 
pitch = PitchFunc(SEApld1d.t);

HeadingFunc = linear_interpolation(nav_unixt[tind], SEAnav1d.Heading[tind], extrapolation_bc=Line());
heading = HeadingFunc(SEApld1d.t);

RollFunc = linear_interpolation(nav_unixt[tind], SEAnav1d.Roll[tind], extrapolation_bc=Line());
roll = RollFunc(SEApld1d.t);

DeclinationFunc = linear_interpolation(nav_unixt[tind], SEAnav1d.Declination[tind], extrapolation_bc=Line());
declination = DeclinationFunc(SEApld1d.t);

DeadReckoningFunc = extrapolate(interpolate((SEAnav1d.t,), Float64.(SEAnav1d.DeadReckoning), Gridded(Constant())), Flat());
deadreckoning = Int64.(DeadReckoningFunc(SEApld1d.t));

seadf = DataFrame(
    time = Dates.unix2datetime.(SEApld1d.t), 
    dive_number = SEApld1d.yo, 
    longitude = SEApld1d.lon,
    latitude = SEApld1d.lat, 
    pressure = SEApld1d.p, 
    nav_resource = SEApld1d.ns,
    declination = declination, 
    pitch = pitch,
    heading = heading,
    roll = roll,
    dead_reckoning = deadreckoning
    );
csvdir = "./"
CSV.write(csvdir * glidername * "_M" * string(missionID) * ".csv", seadf);

=#