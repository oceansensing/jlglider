# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#
# setup directories
rootdir = "/Users/gong/oceansensing Dropbox/C2PO/MARACOOS/electa-20230320-maracoos/";
fromgliderdir = rootdir * "from-glider/"; 
datadir = fromgliderdir * "electa-from-glider-20230405T023652/";
cacdir = fromgliderdir * "cache/";

figoutdir = rootdir * "figures/";
mission = "MARACOOS";
glider = "electa";

# specify valid data time period
t0 = DateTime("2023-03-21");
tN = DateTime("2023-04-21");
trange = datetime2unix.([t0; tN]);

pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 3; # day intervals for plotting
datamode = "realtime"

include("load_slocum_glider.jl")
include("plot_slocum_glider.jl")