# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#
# setup directories

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations
import slocumType: ctdStruct, sciStruct
import slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_ctd_load, glider_presfunc

using Dates
import slocumType: plotSetting, ctdStruct, sciStruct
import load_slocum_glider: load_glider_ctd, load_glider_sci
#import plot_slocum_glider: plot_glider_ctd

datamode = "delayed"; # delayed or realtime
mission = "PASSENGERS 2023";

# specify valid data time period
t0 = DateTime("2023-05-23");
tN = DateTime("2023-06-20");
trange = datetime2unix.([t0; tN]);

pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 1; # day intervals for plotting
ms = 4;
tsms = 4;
pres = (1200, 800)
tspres = (1000, 1000)
ps = plotSetting(pint, iday, ms, tsms, pres, tspres);

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

datamode_sylvia = "delayed"
glidername_sylvia = "sylvia";
rootdir_sylvia = dataroot * "PASSENGERS/2023_glider_data/sylvia-20230608-passengers/";
fromgliderdir_sylvia = rootdir_sylvia * "from-glider/"; 
datadir_sylvia = fromgliderdir_sylvia * datamode * "/" * "sylvia-from-glider-20230612T023321/";
cacdir_sylvia = fromgliderdir_sylvia * "cache/";
figoutdir_sylvia = rootdir_sylvia * "figures/";

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

datadir = datadir_electa;
cacdir = cacdir_electa;
datamode = datamode_electa;
glidername = glidername_electa;

electaCTDraw = load_glider_ctd(datadir_electa, cacdir_electa, trange, datamode_electa, mission, glidername_electa, 1);
#electaCTD = load_glider_ctd(datadir_electa, cacdir_electa, trange, datamode_electa, mission, glidername_electa);
#plot_glider_ctd(electaCTD, figoutdir, ps)

#electaCHLA, electaCDOM, electaBB700, electaBPAR = load_glider_sci(datadir_electa, cacdir_electa, trange, datamode, mission, glidername_electa);
#gliderCHLA, gliderCDOM, gliderBB700, gliderBPAR = electaCHLA, electaCDOM, electaBB700, electaBPAR;

gliderCTD = electaCTDraw;
figoutdir = figoutdir_electa;
include("plot_slocum_glider_ctd.jl")
#include("plot_slocum_glider_bio.jl")


#sylviaCTDraw = load_glider_ctd(datadir_sylvia, cacdir_sylvia, trange, datamode, mission, glidername_sylvia, 1);
#sylviaCHLA, sylviaCDOM, sylviaBB700, sylviaBPAR = load_glider_sci(datadir_sylvia, cacdir_sylvia, trange, datamode, mission, glidername_sylvia);
#gliderCTD = sylviaCTDraw;
#gliderCHLA, gliderCDOM, gliderBB700, gliderBPAR = sylviaCHLA, sylviaCDOM, sylviaBB700, sylviaBPAR;
#figoutdir = figoutdir_sylvia;
#include("plot_slocum_glider_ctd.jl")

#=
nrl641CTD = load_glider_ctd(datadir_nrl641, cacdir_nrl641, trange, datamode, mission, glidername_nrl641);
nrl641CHLA, nrl641CDOM, nrl641BB700, nrl641BPAR = load_glider_sci(datadir_nrl641, cacdir_nrl641, trange, datamode, mission, glidername_nrl641);
gliderCTD = nrl641CTD;
gliderCHLA, gliderCDOM, gliderBB700, gliderBPAR = nrl641CHLA, nrl641CDOM, nrl641BB700, nrl641BPAR;
figoutdir = figoutdir_nrl641;
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
=#