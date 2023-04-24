#using Glob, JLD2, Makie
#import MR_types: MicroRiderRaw
#import MR_io: MR_datasetup, MR_mat2jld2, MR_load_profile
import MR_io: MR_load_profile

project = "NORSE"
mission = "LBE"
profileid = 230;

mrr = MR_load_profile(project, mission, profileid);

