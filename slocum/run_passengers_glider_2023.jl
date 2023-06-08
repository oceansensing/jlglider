# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#
# setup directories
using Dates
import slocumType: plotSetting
import load_slocum_glider: load_glider_ctd
#import plot_slocum_glider: plot_glider_ctd

datamode = "realtime" # delayed or realtime
pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 1; # day intervals for plotting
ms = 8;
tsms = 5;
pres = (1200, 800)
tspres = (1000, 1000)
ps = plotSetting(pint, iday, ms, tsms, pres, tspres);

rootdir = "/Users/gong/oceansensing Dropbox/C2PO/PASSENGERS/2023_glider_data/electa-20230523-passengers/";
fromgliderdir = rootdir * "from-glider/"; 
datadir = fromgliderdir * datamode * "/" * "electa-from-glider-20230608T052727/";
cacdir = fromgliderdir * "cache/";
figoutdir = rootdir * "figures/";
mission = "PASSENGERS 2023";
glider = "electa";

# specify valid data time period
t0 = DateTime("2023-05-23");
tN = DateTime("2023-06-20");
trange = datetime2unix.([t0; tN]);

electaCTD = load_glider_ctd(datadir, cacdir, trange, datamode, mission, glider);
#plot_glider_ctd(electaCTD, figoutdir, ps)

#include("plot_slocum_glider.jl")
