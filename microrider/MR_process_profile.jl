using Glob, JLD2
import MR_types: MicroRiderRaw

datadirJM = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/"
datadirLBE = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/"

datadir = datadirLBE;
pdir = datadir * "mr1000g/";
jld2dir = datadir * "mr1000g_processing/jld2files/";

datafiles = Glob.glob("*.jld2", jld2dir);

i = 1
mrr = MicroRiderRaw[];
data = jldopen(datafiles[i]);
mrr = data["mrprofile"];

