import seaexplorer_functions: load_NAV, load_PLD, seaexplorer_load_mission, seaexplorer_process

mission = 38; # M37 is Jan Mayen, M38 is Lofoten Basin
#include("seaexplorer_load.jl");
jmnav, jmnav1d, jmpld, jmpld1d = seaexplorer_load_mission(37);
lbenav, lbenav1d, lbepld, lbepld1d = seaexplorer_load_mission(38);
jm = seaexplorer_process(jmpld1d);
lbe = seaexplorer_process(lbepld1d);

sea064data = jm;

#include("seaexplorer_plotFast.jl")