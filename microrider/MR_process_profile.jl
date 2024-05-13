# MR_process_profile loads specific MicroRider profiles from NORSE project

include("MR_func.jl")
include("MR_types.jl")

import .MR_types: MicroRiderRaw, MicroRider
#using .MR_func: MR_datasetup, MR_mat2jld2
using .MR_func: MR_load_profile, MR_datasetup

using NaNMath, Statistics, GLMakie, Plots, ColorSchemes, GibbsSeaWater

plotly()

gsw = GibbsSeaWater;
global mr = MR_types.MicroRiderRaw[]; # NOTE!!! if declaring a composite type array as global, then the module name must be included instantiation.

project = "NORSE"
mission = "JM"
year = 2023
#profileid = 1:511;
#profileid = 400:511;
profileid = [5; 154; 427];
#profileid = 331:332;
glider = "SEA064"

reloadflag = 1

#for ii = 1:length(profileid)
for pid in profileid
    if ((@isdefined mr) != true) | (reloadflag == 1)
        display("Loading project " * project * ", mission " * mission * ", profile " * string(pid) * ".")
        mrp = MR_load_profile(project, mission, year, pid); # microrider profile
        push!(mr, mrp);
    end
end

#include("MR_despike_profile.jl")
#include("MR_plots.jl")