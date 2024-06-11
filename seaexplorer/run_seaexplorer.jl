workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

include("seaexplorerFunc.jl")

import .seaexplorerFunc: load_NAV, load_PLD, seaexplorer_load_mission, seaexplorer_process

mission = 48; # M37 is Jan Mayen in 2022, M38 is Lofoten Basin in 2022, M48 is Jan Mayen in 2023
#include("seaexplorer_load.jl");
jm22nav, jm22nav1d, jm22pld, jm22pld1d = seaexplorer_load_mission(37);
lbe22nav, lbe22nav1d, lbe22pld, lbe22pld1d = seaexplorer_load_mission(38);
jm23nav, jm23nav1d, jm23pld, jm23pld1d = seaexplorer_load_mission(48);
jm22 = seaexplorer_process(jm22pld1d);
lbe22 = seaexplorer_process(lbe22pld1d);
jm23 = seaexplorer_process(jm23pld1d);

#sea064pld1d = jm23;
jm = jm23;
lbe = lbe22;

#include("freya_MR_laur_load.jl");
#include("seaexplorer_processing.jl");
#include("seaexplorer_plotFast.jl")
#include("seaexplorer_plotMap.jl")
#include("seaexplorer_plotADCP.jl")
#include("seaexplorer_plotMR.jl")


