# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#
# setup directories

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations

include("slocumType.jl")
include("slocumFunc.jl")
include("slocumLoad.jl")
include("slocumPlot.jl")

import .slocumType: plotSetting, plotStruct, ctdStruct, sciStruct
import .slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_ctd_load, glider_presfunc
import .slocumLoad: load_glider_ctd, load_glider_sci
import .slocumPlot: plot_glider_ctd

datamode = "delayed"; # delayed or realtime
mission = "NORSE-2023";

# specify valid data time period
t0 = DateTime("2023-11-12");
tN = DateTime("2023-11-27");
trange = datetime2unix.([t0; tN]);

lonrange = [-9.5 -6.0];
latrange = [70.5 71.5]; 

pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 1; # day intervals for plotting
ms = 4; # marker size
tsms = 4; # time series marker size
pres = (1200, 800); # plot resolution
tspres = (1000, 1000); # time series plot resolution
fs = 24; # font size
ps = plotSetting(pint, iday, ms, tsms, pres, tspres, fs);

#dataroot = "/mnt/c/Users/C2PO/oceansensing Dropbox/C2PO/";
dataroot = "/Users/gong/oceansensing Dropbox/C2PO/";

datamode_electa = "delayed"
glidername_electa = "electa";
rootdir_electa = dataroot * "glider/gliderData/electa-20231127-norse-post/";

fromgliderdir_electa = rootdir_electa; 
if datamode_electa == "delayed"
    datadir_electa = fromgliderdir_electa * "delayed/";
elseif datamode_electa == "realtime"
    datadir_electa = fromgliderdir_electa * "realtime/";
end
cacdir_electa = fromgliderdir_electa * "cache/";
figoutdir_electa = rootdir_electa * "figures/";
loadmode_electa = "lowercase";
ctemprange = (-1, 5);
condrange = (2.65, 3.4);
saltArange = (33.5, 35.2);
sigma0range = (26.4, 28.4);
sndspdrange = (1450, 1470);
spice0range = (-1, 1);
temprange = ctemprange;
saltrange = saltArange;
pst_electa = plotStruct(figoutdir_electa, mission, glidername_electa, temprange[1], temprange[2], condrange[1], condrange[2], saltrange[1], saltrange[2], sigma0range[1], sigma0range[2], spice0range[1], spice0range[2], sndspdrange[1], sndspdrange[2]);

datadir = datadir_electa;
cacdir = cacdir_electa;
datamode = datamode_electa;
glidername = glidername_electa;
loadmode = loadmode_electa;

#electaCTDraw = load_glider_ctd(datadir_electa, cacdir_electa, trange, lonrange, latrange, datamode_electa, mission, glidername_electa, loadmode_electa);
plot_glider_ctd(electaCTDraw, ps, pst_electa);

#include("slocumLoadTest.jl")

pst = pst_electa;
glider1 = electaCTDraw;
glider2 = glider1;
#plot_glider_map(electaCTDraw, sylviaCTDraw, ps, pst);
include("slocumPlotMapTest.jl");