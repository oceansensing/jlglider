using Glob, JLD2
import MR_types: MicroRiderRaw
import MR_io: MR_mat2jld2, MR_loadjld2, MR_datasetup

jld2dir, matdir, pdir = MR_datasetup("NORSE", "LBE");
datafiles = Glob.glob("*.jld2", jld2dir);

i = 1
mrr = MR_loadjld2(datafiles[1]);