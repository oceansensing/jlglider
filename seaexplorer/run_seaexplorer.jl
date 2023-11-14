workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

import seaexplorer_functions: load_NAV, load_PLD, seaexplorer_load_mission, seaexplorer_process

mission = 48; # M37 is Jan Mayen in 2022, M38 is Lofoten Basin in 2022, M48 is Jan Mayen in 2023
#include("seaexplorer_load.jl");
jmnav, jmnav1d, jmpld, jmpld1d = seaexplorer_load_mission(mission);
lbenav, lbenav1d, lbepld, lbepld1d = seaexplorer_load_mission(mission);
jm = seaexplorer_process(jmpld1d);
lbe = seaexplorer_process(lbepld1d);
#sea064data = jm;

include("freya_MR_laur_load.jl");

#include("seaexplorer_processing.jl");
#include("seaexplorer_plotFast.jl")

