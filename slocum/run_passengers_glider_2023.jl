# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#
# setup directories

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations
import slocumType: plotSetting, plotStruct, ctdStruct, sciStruct
import slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_ctd_load, glider_presfunc
import slocumLoad: load_glider_ctd, load_glider_sci
import slocumPlot: plot_glider_ctd

datamode = "delayed"; # delayed or realtime
mission = "NESMA-PASSENGERS-2023";

# specify valid data time period
t0 = DateTime("2023-05-23");
tN = DateTime("2023-06-20");
trange = datetime2unix.([t0; tN]);

lonrange = [-65.0 -55.0];
latrange = [30.0 50.0]; 

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
rootdir_electa = dataroot * "PASSENGERS/2023_glider_data/electa-20230523-passengers/";
fromgliderdir_electa = rootdir_electa * "from-glider/"; 
if datamode_electa == "delayed"
    datadir_electa = fromgliderdir_electa * datamode_electa * "/";
elseif datamode_electa == "realtime"
    datadir_electa = fromgliderdir_electa * datamode_electa * "/electa/from-glider/";
end
cacdir_electa = fromgliderdir_electa * "cache/";
figoutdir_electa = rootdir_electa * "figures/";
loadmode_electa = "lowercase";
temprange = (19, 25);
condrange = (4.9, 5.5);
saltrange = (36.3, 37.1);
sigma0range = (24.5, 26.4);
sndspdrange = (1522, 1536);
spice0range = (4.2, 6.0);
pst_electa = plotStruct(figoutdir_electa, mission, glidername_electa, temprange[1], temprange[2], condrange[1], condrange[2], saltrange[1], saltrange[2], sigma0range[1], sigma0range[2], spice0range[1], spice0range[2], sndspdrange[1], sndspdrange[2]);


datamode_sylvia = "delayed"
glidername_sylvia = "sylvia";
rootdir_sylvia = dataroot * "PASSENGERS/2023_glider_data/sylvia-20230608-passengers/";
fromgliderdir_sylvia = rootdir_sylvia * "from-glider/";
if datamode_sylvia == "delayed"
    datadir_sylvia = fromgliderdir_sylvia * datamode_sylvia * "/";
elseif datamode_sylvia == "realtime"
    datadir_sylvia = fromgliderdir_sylvia * datamode_sylvia * "/" * glidername_sylvia * "/from-glider/";
end
#datadir_sylvia = fromgliderdir_sylvia * datamode * "/" * "sylvia-from-glider-20230612T023321/";
cacdir_sylvia = fromgliderdir_sylvia * "cache/";
figoutdir_sylvia = rootdir_sylvia * "figures/";
loadmode_sylvia = "uppercase";
temprange = (19, 25);
condrange = (4.9, 5.5);
saltrange = (36.3, 37.1);
sigma0range = (24.5, 26.4);
sndspdrange = (1522, 1536);
spice0range = (4.2, 6.0);
pst_sylvia = plotStruct(figoutdir_sylvia, mission, glidername_sylvia, temprange[1], temprange[2], condrange[1], condrange[2], saltrange[1], saltrange[2], sigma0range[1], sigma0range[2], spice0range[1], spice0range[2], sndspdrange[1], sndspdrange[2]);

glidername_nrl641 = "nrl641";
rootdir_nrl641 = dataroot * "PASSENGERS/2023_glider_data/nrl641-20230523-passengers/";
fromgliderdir_nrl641 = rootdir_nrl641 * "from-glider/"; 
datadir_nrl641 = fromgliderdir_nrl641 * datamode * "/" * "f641sg17-from-glider-20230610T113712/";
cacdir_nrl641 = fromgliderdir_nrl641 * "cache/";
figoutdir_nrl641 = rootdir_nrl641 * "figures/";

glidername_ru30 = "ru30";
rootdir_ru30 = dataroot * "PASSENGERS/2023_glider_data/ru30-20230525-passengers/";
fromgliderdir_ru30 = rootdir_ru30 * "from-glider/"; 
datadir_ru30 = fromgliderdir_ru30 * datamode * "/" * "ru30/from-glider/";
cacdir_ru30 = fromgliderdir_ru30 * "cache/";
figoutdir_ru30 = rootdir_ru30 * "figures/";

glidername_ru36 = "ru36";
rootdir_ru36 = dataroot * "PASSENGERS/2023_glider_data/ru36-20230525-passengers/";
fromgliderdir_ru36 = rootdir_ru36 * "from-glider/"; 
datadir_ru36 = fromgliderdir_ru36 * datamode * "/" * "ru36/from-glider/";
cacdir_ru36 = fromgliderdir_ru36 * "cache/";
figoutdir_ru36 = rootdir_ru36 * "figures/";

datadir = datadir_sylvia;
cacdir = cacdir_sylvia;
datamode = datamode_sylvia;
glidername = glidername_sylvia;
loadmode = loadmode_sylvia;

electaCTDraw = load_glider_ctd(datadir_electa, cacdir_electa, trange, lonrange, latrange, datamode_electa, mission, glidername_electa, loadmode_electa);
#plot_glider_ctd(electaCTDraw, ps, pst_electa);

sylviaCTDraw = load_glider_ctd(datadir_sylvia, cacdir_sylvia, trange, lonrange, latrange, datamode_sylvia, mission, glidername_sylvia, loadmode_sylvia);
#include("slocumLoadTest.jl")
#plot_glider_ctd(sylviaCTDraw, ps, pst_sylvia);

pst = pst_electa;
glider1 = electaCTDraw;
glider2 = sylviaCTDraw;
#plot_glider_map(electaCTDraw, sylviaCTDraw, ps, pst);
include("slocumPlotMapTest.jl");