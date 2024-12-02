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
import .gliderPlot: plot_glider_ctd, plotGliderCTD, plotGliderMap
using JLD2

reloadflag = false

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdirpath = "/Users/gong/GitHub/jlglider/seaexplorer/seaexplorer_yaml_NORSE/";

#SEAnav, SEAnav1d, SEApld, SEApld1d = seaexplorer_load_mission(missionYAMLdirpath * "sea094-20240709-passengers.yaml");

if @isdefined(gliderCTDarray) == false
    display("Loading data...")
    #global gliderCTDarray = SeaExplorerCTD[];
    if reloadflag == true
        gliderCTDarray = seaexplorerYAMLload(missionYAMLdirpath);
        jldsave(gliderdatadir * "NORSE_seaexplorerCTDdata.jld2"; gliderCTDarray);
        display("Done reloading data.")
    else
        gliderCTDarray = load(gliderdatadir * "NORSE_seaexplorerCTDdata.jld2")["gliderCTDarray"];
        display("Done loading data.")
    end
end

global ps = Glider.gliderPlotType.plotSetting[];
global pst = Glider.gliderPlotType.plotStruct[];
for i = 1:length(gliderCTDarray)
    pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
    iday = 1; # day intervals for plotting
    ms = 10; # marker size
    tsms = 6; # time series marker size
    pres = (1600, 500); # plot resolution
    tspres = (1000, 1000); # time series plot resolution
    fs = 32; # font size
    global ps = push!(ps, Glider.gliderPlotType.plotSetting(pint, iday, ms, tsms, pres, tspres, fs));

    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/NORSE/";
    project = gliderCTDarray[i].project;
    glidername = gliderCTDarray[i].glidername;
    missionID = gliderCTDarray[i].missionID;
    if missionID == 38
        lonmin, lonmax = -10, -5.5;    
        latmin, latmax = 70.2, 71.5;
        zlo, zhi = -1000, 0;
        tempmin, tempmax = -2, 9;
        condmin, condmax = 25, 37;
        saltmin, saltmax = 35, 35.4;
        sigma0min, sigma0max = 26.8, 28.2;
        spice0min, spice0max = -1.25, 1.25;
        sndspdmin, sndspdmax = 1450, 1485;
    else
        lonmin, lonmax = -10, -5.5;    
        latmin, latmax = 70.2, 71.5;
        zlo, zhi = -1000, 0;
        tempmin, tempmax = -1, 6.5;
        condmin, condmax = 27, 34;
        saltmin, saltmax = 33.7, 35.2;
        sigma0min, sigma0max = 26.8, 28.2;
        spice0min, spice0max = -1.0, 0.55;
        sndspdmin, sndspdmax = 1450, 1475;
    end
    global pst = push!(pst, Glider.gliderPlotType.plotStruct(figoutdir, project, glidername, lonmin, lonmax, latmin, latmax, zlo, zhi, tempmin, tempmax, condmin, condmax, saltmin, saltmax, sigma0min, sigma0max, spice0min, spice0max, sndspdmin, sndspdmax));
end

plotGliderCTD(gliderCTDarray, ps, pst)
plotGliderMap(gliderCTDarray, pst, pzrange=[-20,-10], varname="spice0", logzflag=0);

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


