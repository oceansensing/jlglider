# this script loads Slocum glider data using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
# gong@vims.edu 2024-09-05: added a function to load glider data from yaml metadata files for a general mission
# gong@vims.edu 2024-11-11: major refactoring to unify plotting for seaexplorer and slocum glider data

workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

workdir = "/Users/gong/GitHub/jlglider/slocum"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

using JLD2, Glider

include("/Users/gong/GitHub/jlglider/slocum/slocumLoad.jl")
include("/Users/gong/GitHub/ocean_julia/C2PO.jl")
include("/Users/gong/GitHub/jlglider/seaexplorer/gliderPlot.jl")

import .slocumLoad: load_glider_ctd, load_glider_sci, glider_ctd_qc, slocumYAMLload
import .gliderPlot: plotGliderCTD, plotGliderMap

reloadflag = false

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdirpath = "/Users/gong/GitHub/jlglider/slocum/mission_yaml_MAB/sylvia-20240926-hurricane.yaml";

if @isdefined(gliderCTDarray) == false
    if reloadflag == true
        gliderCTDarray = slocumYAMLload(missionYAMLdirpath);
        jldsave(gliderdatadir * "MARACOOS_slocumCTDdata_20240926.jld2"; gliderCTDarray);
        display("Done reloading data.")
    else
        gliderCTDarray = load(gliderdatadir * "MARACOOS_slocumCTDdata_20240926.jld2")["gliderCTDarray"];
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
    pres = (1600, 800); # plot resolution
    tspres = (1000, 1000); # time series plot resolution
    fs = 32; # font size
    global ps = push!(ps, Glider.gliderPlotType.plotSetting(pint, iday, ms, tsms, pres, tspres, fs));

    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/MAB/";
    project = gliderCTDarray[i].project;
    glidername = gliderCTDarray[i].glidername;
    latmin, latmax = 35, 42;
    lonmin, lonmax = -77.5, -69;    
    zlo, zhi = -100, 0;
    tempmin, tempmax = 9.0, 23.1;
    condmin, condmax = 3.5, 5.2;
    saltmin, saltmax = 29.0, 35.5;
    sigma0min, sigma0max = 20.9, 26.5;
    spice0min, spice0max = -1, 4;
    sndspdmin, sndspdmax = 1490, 1530;
    global pst = push!(pst, Glider.gliderPlotType.plotStruct(figoutdir, project, glidername, lonmin, lonmax, latmin, latmax, zlo, zhi, tempmin, tempmax, condmin, condmax, saltmin, saltmax, sigma0min, sigma0max, spice0min, spice0max, sndspdmin, sndspdmax));
end

plotGliderCTD(gliderCTDarray, ps, pst)
plotGliderMap(gliderCTDarray, pst, pzrange=[-40,-30], varname="saltA", logzflag=0);
