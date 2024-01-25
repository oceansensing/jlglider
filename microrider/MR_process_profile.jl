# MR_process_profile loads specific MicroRider profiles from NORSE project

include("MR_io.jl")
include("MR_types.jl")

import .MR_types: MicroRiderRaw
#using .MR_io: MR_datasetup, MR_mat2jld2
using .MR_io: MR_load_profile

using NaNMath, GLMakie, ColorSchemes, GibbsSeaWater

gsw = GibbsSeaWater;
global mrp = MicroRiderRaw[];

project = "NORSE"
mission = "JM"
year = 2023
profileid = 101
#profileid = 150
glider = "SEA064"

reloadflag = 1

if ((@isdefined mrr) != true) | (reloadflag == 1)
    display("Loading project " * project * ", mission " * mission * ", profile " * string(profileid) * ".")
    mrp = MR_load_profile(project, mission, year, profileid); # microrider profile
    display("Loaded project " * project * ", mission " * mission * ", profile " * string(profileid) * ".");
end

mrpz = gsw.gsw_z_from_p.(mrp.P_fast, 71.0, 0.0, 0.0);

include("MR_plots.jl")