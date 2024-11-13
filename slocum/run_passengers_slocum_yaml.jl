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

reloadflag = true

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdir = "/Users/gong/GitHub/jlglider/slocum/slocum_yaml_PASSENGERS/";

if @isdefined(gliderCTDarray) == false
    if reloadflag == true
        gliderCTDarray = slocumYAMLload(missionYAMLdir);
        jldsave(gliderdatadir * "PASSENGERS_slocumCTDdata.jld2"; gliderCTDarray);
        display("Done reloading data.")
    else
        gliderCTDarray = load(gliderdatadir * "PASSENGERS_slocumCTDdata.jld2")["gliderCTDarray"];
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

    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/NESMA-PASSENGERS/";
    project = gliderCTDarray[i].project;
    glidername = gliderCTDarray[i].glidername;
    latmin, latmax = 37, 40;
    lonmin, lonmax = -65.2, -59.8;    
    tempmin, tempmax = 4.0, 32.0;
    condmin, condmax = 3, 6.5;
    saltmin, saltmax = 31.0, 37.25;
    sigma0min, sigma0max = 20.0, 30.0;
    spice0min, spice0max = -1.5, 7.5;
    sndspdmin, sndspdmax = 1480, 1550;
    global pst = push!(pst, Glider.gliderPlotType.plotStruct(figoutdir, project, glidername, lonmin, lonmax, latmin, latmax, tempmin, tempmax, condmin, condmax, saltmin, saltmax, sigma0min, sigma0max, spice0min, spice0max, sndspdmin, sndspdmax));
end

plotGliderCTD(gliderCTDarray, ps, pst)
plotGliderMap(gliderCTDarray, pst, pzrange=[-40,-30], varname="saltA", logzflag=0);
