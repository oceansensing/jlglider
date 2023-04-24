using Glob, JLD2
import MR_types: MicroRiderRaw
import MR_io: MR_mat2jld2, MR_loadjld2, MR_datasetup, MR_load_profile

project = "NORSE"
mission = "LBE"

mrr = MR_load_profile(project, mission, 15);