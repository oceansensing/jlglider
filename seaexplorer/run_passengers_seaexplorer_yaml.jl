workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

# loading the Glider module with all relevant types
using Glider

# import all the functions used
include("/Users/gong/GitHub/jlglider/seaexplorer/seaexplorerFunc.jl")
include("/Users/gong/GitHub/ocean_julia/C2PO.jl")
include("/Users/gong/GitHub/jlglider/seaexplorer/gliderPlot.jl")
import .seaexplorerFunc: seaexplorerYAMLload, seaexplorer_load_mission 
import .gliderPlot: plot_glider_ctd, plotGliderCTD
using JLD2

reloadflag = false

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdirpath = "/Users/gong/GitHub/jlglider/seaexplorer/mission_yaml_PASSENGERS/";

#SEAnav, SEAnav1d, SEApld, SEApld1d = seaexplorer_load_mission(missionYAMLdirpath * "sea094-20240709-passengers.yaml");

if @isdefined(gliderCTDarray) == false
    display("Loading data...")
    #global gliderCTDarray = SeaExplorerCTD[];
    if reloadflag == true
        gliderCTDarray = seaexplorerYAMLload(missionYAMLdirpath);
        jldsave(gliderdatadir * "PASSENGERS_seaexplorerCTDdata.jld2"; gliderCTDarray);
        display("Done reloading data.")
    else
        gliderCTDarray = load(gliderdatadir * "PASSENGERS_seaexplorerCTDdata.jld2")["gliderCTDarray"];
        display("Done loading data.")
    end
end

global ps = Glider.gliderPlotType.plotSetting[];
global pst = Glider.gliderPlotType.plotStruct[];
for i = 1:length(gliderCTDarray)
    pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
    iday = 1; # day intervals for plotting
    ms = 6; # marker size
    tsms = 6; # time series marker size
    pres = (1600, 800); # plot resolution
    tspres = (1000, 1000); # time series plot resolution
    fs = 42; # font size
    global ps = push!(ps, Glider.gliderPlotType.plotSetting(pint, iday, ms, tsms, pres, tspres, fs));

    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/";
    project = gliderCTDarray[i].project;
    glidername = gliderCTDarray[i].glidername;
    tempmin = 5.0;
    tempmax = 30.0;
    condmin = 30;
    condmax = 50;
    saltmin = 34.5;
    saltmax = 37.5;
    sigma0min = 20.0;
    sigma0max = 30.0;
    spice0min = -5.0;
    spice0max = 5.0;
    sndspdmin = 1450;
    sndspdmax = 1530;
    global pst = push!(pst, Glider.gliderPlotType.plotStruct(figoutdir, project, glidername, tempmin, tempmax, condmin, condmax, saltmin, saltmax, sigma0min, sigma0max, spice0min, spice0max, sndspdmin, sndspdmax));
end

plotGliderCTD(gliderCTDarray, ps, pst)


#plotSeaExplorerCTD(gliderCTDarray)

#missionYAML = "sea064-20240720-nesma.yaml";
#sea064nav, sea064nav1d, sea064pld, sea064pld1d = seaexplorer_load_mission(missionYAMLdir * missionYAML);
#sea064nesma0720 = seaexplorer_process(sea064pld1d);

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


