# MR_process_profile loads specific MicroRider profiles from NORSE project

include("MR_io.jl")
include("MR_types.jl")

import .MR_types: MicroRiderRaw, MicroRider
#using .MR_io: MR_datasetup, MR_mat2jld2
using .MR_io: MR_load_profile

using NaNMath, GLMakie, ColorSchemes, GibbsSeaWater

gsw = GibbsSeaWater;
global norse23mr = MicroRider[];

project = "NORSE"
mission = "JM"
year = 2023
#profileid = 1:511;
profileid = 150;
#profileid = 150
glider = "SEA064"

reloadflag = 1

for ii = 1:length(profileid)
    if ((@isdefined norse23mr) != true) | (reloadflag == 1)
        display("Loading project " * project * ", mission " * mission * ", profile " * string(profileid[ii]) * ".")
        mrp = MR_load_profile(project, mission, year, profileid[ii]); # microrider profile
        mrpz = gsw.gsw_z_from_p.(mrp.P_fast, 71.0, 0.0, 0.0); 
        global norse23mr = push!(norse23mr, MicroRider(mrp, mrpz));
        #display("Loaded project " * project * ", mission " * mission * ", profile " * string(profileid[ii]) * ".");
    end
end

#mrpz = gsw.gsw_z_from_p.(mrp.P_fast, 71.0, 0.0, 0.0);

include("MR_despike_profile.jl")
#include("MR_plots.jl")