using NaNMath, GLMakie, ColorSchemes

#import MR_types: MicroRiderRaw
#import MR_io: MR_datasetup, MR_mat2jld2, MR_load_profile
import MR_io: MR_load_profile

project = "NORSE"
mission = "LBE"
profileid = 2
glider = "SEA064"

reloadflag = 1

if ((@isdefined mrr) != true) | (reloadflag == 1)
    display("Loading project " * project * ", mission " * mission * ", profile " * string(profileid) * ".")
    mrr = MR_load_profile(project, mission, profileid);
end

#include("MR_plots.jl")