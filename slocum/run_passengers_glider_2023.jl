# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#
# setup directories
using Dates
import slocumType: plotSetting, ctdStruct, sciStruct
import load_slocum_glider: load_glider_ctd, load_glider_sci
#import plot_slocum_glider: plot_glider_ctd

datamode = "realtime" # delayed or realtime
pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 1; # day intervals for plotting
ms = 8;
tsms = 5;
pres = (1200, 800)
tspres = (1000, 1000)
ps = plotSetting(pint, iday, ms, tsms, pres, tspres);

mission = "PASSENGERS 2023";

glidername_electa = "electa";
rootdir_electa = "/Users/gong/oceansensing Dropbox/C2PO/PASSENGERS/2023_glider_data/electa-20230523-passengers/";
fromgliderdir_electa = rootdir_electa * "from-glider/"; 
datadir_electa = fromgliderdir_electa * datamode * "/" * "electa-from-glider-20230609T181801/";
cacdir_electa = fromgliderdir_electa * "cache/";
figoutdir_electa = rootdir_electa * "figures/";

glidername_nrl641 = "nrl641";
rootdir_nrl641 = "/Users/gong/oceansensing Dropbox/C2PO/PASSENGERS/2023_glider_data/nrl641-20230523-passengers/";
fromgliderdir_nrl641 = rootdir_nrl641 * "from-glider/"; 
datadir_nrl641 = fromgliderdir_nrl641 * datamode * "/" * "f641sg17-from-glider-20230610T040435/";
cacdir_nrl641 = fromgliderdir_nrl641 * "cache/";
figoutdir_nrl641 = rootdir_nrl641 * "figures/";

glidername_ru30 = "ru30";
rootdir_ru30 = "/Users/gong/oceansensing Dropbox/C2PO/PASSENGERS/2023_glider_data/ru30-20230525-passengers/";
fromgliderdir_ru30 = rootdir_ru30 * "from-glider/"; 
datadir_ru30 = fromgliderdir_ru30 * datamode * "/" * "ru30-from-glider-20230609T185136/";
cacdir_ru30 = fromgliderdir_ru30 * "cache/";
figoutdir_ru30 = rootdir_ru30 * "figures/";

glidername_ru36 = "ru36";
rootdir_ru36 = "/Users/gong/oceansensing Dropbox/C2PO/PASSENGERS/2023_glider_data/ru36-20230525-passengers/";
fromgliderdir_ru36 = rootdir_ru36 * "from-glider/"; 
datadir_ru36 = fromgliderdir_ru36 * datamode * "/" * "ru36-from-glider-20230609T182417/";
cacdir_ru36 = fromgliderdir_ru36 * "cache/";
figoutdir_ru36 = rootdir_ru36 * "figures/";

# specify valid data time period
t0 = DateTime("2023-05-23");
tN = DateTime("2023-06-20");
trange = datetime2unix.([t0; tN]);

electaCTD = load_glider_ctd(datadir_electa, cacdir_electa, trange, datamode, mission, glidername_electa);
electaCHLA, electaCDOM, electaBB700, electaBPAR = load_glider_sci(datadir_electa, cacdir_electa, trange, datamode, mission, glidername_electa);
#plot_glider_ctd(electaCTD, figoutdir, ps)
gliderCTD = electaCTD;
gliderCHLA, gliderCDOM, gliderBB700, gliderBPAR = electaCHLA, electaCDOM, electaBB700, electaBPAR;
figoutdir = figoutdir_electa;
include("plot_slocum_glider.jl")

ru30CTD = load_glider_ctd(datadir_ru30, cacdir_ru30, trange, datamode, mission, glidername_ru30);
ru30CHLA, ru30CDOM, ru30BB700, ru30BPAR = load_glider_sci(datadir_ru30, cacdir_ru30, trange, datamode, mission, glidername_ru30);
gliderCTD = ru30CTD;
gliderCHLA, gliderCDOM, gliderBB700, gliderBPAR = ru30CHLA, ru30CDOM, ru30BB700, ru30BPAR;
figoutdir = figoutdir_ru30;
include("plot_slocum_glider.jl")

ru36CTD = load_glider_ctd(datadir_ru36, cacdir_ru36, trange, datamode, mission, glidername_ru36);
ru36CHLA, ru36CDOM, ru36BB700, ru36BPAR = load_glider_sci(datadir_ru36, cacdir_ru36, trange, datamode, mission, glidername_ru36);
gliderCTD = ru36CTD;
gliderCHLA, gliderCDOM, gliderBB700, gliderBPAR = ru36CHLA, ru36CDOM, ru36BB700, ru36BPAR;
figoutdir = figoutdir_ru36;
include("plot_slocum_glider.jl")
