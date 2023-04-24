using Glob, JLD2
#import MR_types: MicroRiderRaw
#import MR_io: MR_datasetup, MR_mat2jld2, MR_load_profile
import MR_io: MR_load_profile

project = "NORSE"
mission = "LBE"

mrr = MR_load_profile(project, mission, 15);