# this script loads Slocum glider data using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
# gong@vims.edu 2024-09-05: added a function to load glider data from yaml metadata files for a general mission
#

using PyCall
using Glob, YAML, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations, YAML, JLD2

include("slocumType.jl")
include("slocumFunc.jl")
include("slocumLoad.jl")
include("slocumPlot.jl")

using .slocumType: plotSetting, plotStruct, ctdStruct, sciStruct
#using Main.slocumLoad.slocumType: ctdStruct
import .slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_presfunc
import .slocumLoad: load_glider_ctd, load_glider_sci, glider_ctd_qc, slocumYAMLload
import .slocumPlot: plot_glider_ctd

reloadflag = false

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdir = "/Users/gong/GitHub/jlglider/slocum/mission_yaml/";

if @isdefined(gliderCTDarray) == false
    if reloadflag == true
        gliderCTDarray = slocumYAMLload(missionYAMLdir);
        jldsave(gliderdatadir * "slocumCTDdata.jld2"; gliderCTDarray);
        display("Done reloading data.")
    else
        gliderCTDarray = load(gliderdatadir * "slocumCTDdata.jld2")["gliderCTDarray"];
        display("Done loading data.")
    end
end


for i = 1:length(gliderCTDarray)
    #i = 6
    gliderCTDraw = gliderCTDarray[i];
    lonrange = [NaNMath.minimum(gliderCTDraw.lon) NaNMath.maximum(gliderCTDraw.lon)];
    latrange = [NaNMath.minimum(gliderCTDraw.lat) NaNMath.maximum(gliderCTDraw.lat)]; 

    pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
    iday = 1; # day intervals for plotting
    ms = 6; # marker size
    tsms = 6; # time series marker size
    pres = (1200, 800); # plot resolution
    tspres = (1000, 1000); # time series plot resolution
    fs = 24; # font size
    ps = plotSetting(pint, iday, ms, tsms, pres, tspres, fs);

    ctempstd = std(gliderCTDraw.ctemp);
    condstd = std(gliderCTDraw.cond);
    saltAstd = std(gliderCTDraw.saltA);
    sigma0std = std(gliderCTDraw.sigma0);
    sndspdstd = std(gliderCTDraw.sndspd);
    spice0std = std(gliderCTDraw.spice0);

    nsig = 2.5

    #figoutdir = "/Users/gong/GitHub/jlglider/slocum/figures/";
    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/";
    ctemprange = (NaNMath.mean(gliderCTDraw.ctemp) .- nsig*ctempstd, NaNMath.mean(gliderCTDraw.ctemp) .+ nsig*ctempstd);
    condrange = (NaNMath.mean(gliderCTDraw.cond) .- nsig*condstd, NaNMath.mean(gliderCTDraw.cond) .+ nsig*condstd);
    saltArange = (NaNMath.mean(gliderCTDraw.saltA) .- nsig*saltAstd, NaNMath.mean(gliderCTDraw.saltA) .+ nsig*saltAstd);
    sigma0range = (NaNMath.mean(gliderCTDraw.sigma0) .- nsig*sigma0std, NaNMath.mean(gliderCTDraw.sigma0) .+ nsig*sigma0std);
    sndspdrange = (NaNMath.mean(gliderCTDraw.sndspd) .- nsig*sndspdstd, NaNMath.mean(gliderCTDraw.sndspd) .+ nsig*sndspdstd);
    spice0range = (NaNMath.mean(gliderCTDraw.spice0) .- nsig*spice0std, NaNMath.mean(gliderCTDraw.spice0) .+ nsig*spice0std);
    temprange = ctemprange;
    saltrange = saltArange;
    pst = plotStruct(figoutdir, gliderCTDraw.mission, gliderCTDraw.glidername, temprange[1], temprange[2], condrange[1], condrange[2], saltrange[1], saltrange[2], sigma0range[1], sigma0range[2], spice0range[1], spice0range[2], sndspdrange[1], sndspdrange[2]);

    plot_glider_ctd(gliderCTDraw, ps, pst);
end

#include("slocumLoadTest.jl")

#=pst = pst_electa;
glider1 = electaCTDraw;
glider2 = glider1;
#plot_glider_map(electaCTDraw, sylviaCTDraw, ps, pst);
include("slocumPlotMapTest.jl");
=#