workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

include("seaexplorerFunc.jl")

import .seaexplorerFunc: load_NAV, load_PLD, seaexplorer_load_mission, seaexplorer_process

reloadflag = true

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdir = "/Users/gong/GitHub/jlglider/seaexplorer/mission_yaml_NESMA/";

if @isdefined(gliderCTDarray) == false
    if reloadingflag == true
        gliderCTDarray = seaexplorerYAMLload(missionYAMLdir);
        jldsave(gliderdatadir * "PASSENGERS_seaexplorerCTDdata.jld2"; gliderCTDarray);
        display("Done reloading data.")
    else
        gliderCTDarray = load(gliderdatadir * "PASSENGERS_seaexplorerCTDdata.jld2")["gliderCTDarray"];
        display("Done loading data.")
    end
end

#plotSeaExplorerCTD(gliderCTDarray)

missionYAML = "sea064-20240720-nesma.yml";
sea064nav, sea064nav1d, sea064pld, sea064pld1d = seaexplorer_load_mission(missionYAML);
sea064nesma0720 = seaexplorer_process(sea064pld1d);

#missionYAML = "sea094-20240709-nesma.yml";
#sea094nav, sea094nav1d, sea094pld, sea094pld1d = seaexplorer_load_mission(missionYAML);
#sea094nesma0709 = seaexplorer_process(sea064pld1d);


# SEA064: M37 is Jan Mayen in 2022, M38 is Lofoten Basin in 2022, M48 is Jan Mayen in 2023, M58 is NESMA 2024
#sea064nav, sea064nav1d, sea064pld, sea064pld1d = seaexplorer_load_mission("sea064", 58)
# SEA094: M41 is NESMA 2024
#sea094nav, sea094nav1d, sea094pld, sea094pld1d = seaexplorer_load_mission("sea094", 41);
#sea094 = seaexplorer_process(sea094pld1d);

#include("freya_MR_laur_load.jl");
#include("seaexplorer_processing.jl");
#include("seaexplorer_plotFast.jl")
#include("seaexplorer_plotMap.jl")
#include("seaexplorer_plotADCP.jl")
#include("seaexplorer_plotMR.jl")


