using JLD2, Glob
import MR_types: MicroRiderRaw

datadirJM = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/"
datadirLBE = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/"

datadir = datadirLBE;
pdir = datadir * "mr1000g/";
matdir = datadir * "mr1000g_processing/matfiles/";
jld2dir = datadir * "mr1000g_processing/jld2files/";

mrr = MicroRiderRaw[];

if isdir(jld2dir)
    jldfiles = Glob.glob("*.jld2", jld2dir);
    if !isempty(jldfiles)
        i = 1
        display(i)
        data = jldopen(jldfiles[i], "r");
        mrr = data["mrprofile"];
    end #if
end #if
