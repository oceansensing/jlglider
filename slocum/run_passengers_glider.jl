# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#
# setup directories
rootdir = "/Users/gong/oceansensing Dropbox/C2PO/PASSENGERS/2022_glider_data/electa-20221105-passengers/";
fromgliderdir = rootdir * "from-glider/"; 
datadir = fromgliderdir;
cacdir = fromgliderdir * "cache/";

figoutdir = rootdir * "figures/";
mission = "PASSENGERS";
glider = "electa";

# specify valid data time period
t0 = DateTime("2022-10-01");
tN = DateTime("2022-12-01");
trange = datetime2unix.([t0; tN]);

datamode = "delayed"

include("load_slocum_glider.jl")
include("plot_slocum_glider.jl")