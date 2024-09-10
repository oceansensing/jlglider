# this script loads Slocum glider data using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
# gong@vims.edu 2024-09-05: added a function to load glider data from yaml metadata files for a general mission
#

using PyCall
using Glob, YAML, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations, YAML, JLD2, NCDatasets, GLMakie

include("slocumType.jl")
include("slocumFunc.jl")
include("slocumLoad.jl")
include("slocumPlot.jl")

using .slocumType: plotSetting, plotStruct, ctdStruct, sciStruct
#using Main.slocumLoad.slocumType: ctdStruct
import .slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_presfunc
import .slocumLoad: load_glider_ctd, load_glider_sci, glider_ctd_qc, slocumYAMLload
import .slocumPlot: plot_glider_ctd, plotSlocumCTD

reloadflag = false

gliderdatadir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/"; 
missionYAMLdir = "/Users/gong/GitHub/jlglider/slocum/mission_yaml/";

if @isdefined(gliderCTDarray)
    if reloadflag == true
        gliderCTDarray = slocumYAMLload(missionYAMLdir);
        jldsave(gliderdatadir * "slocumCTDdata.jld2"; gliderCTDarray);
        display("Done reloading data.")
    else
        gliderCTDarray = load(gliderdatadir * "slocumCTDdata.jld2")["gliderCTDarray"];
        display("Done loading data.")
    end
end

#plotSlocumCTD(gliderCTDarray)

bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
bathyds = Dataset(bathypath,"r");

lon = bathyds["lon"][:];
lat = bathyds["lat"][:];

latmin, latmax = 35, 42;
lonmin, lonmax = -77.5, -69;

# approximate x-axis scaling to make it look "normal"
dlat = latmax - latmin;
dlon = lonmax - lonmin;
lat0 = NaNMath.mean([38.0]);
xfac = sind(90-lat0);
yres = 2000;
pres = (ceil(yres * (xfac/(dlat/dlon))), yres);

# extract indices from the bathymetric data file
latind = findall(latmin-0.1 .<= lat .<= latmax+0.1);
lonind = findall(lonmin-0.1 .<= lon .<= lonmax+0.1);

z = Float64.(bathyds["z"][lonind, latind]); # recasting as Float64 to fix a StackOverFlow error seen in GLMakie 0.6.0.
x = lon[lonind];
y = lat[latind];

pzind = findall(z .> 1);
nzind = findall(z .< -1);
zzind = findall(-1 .<= z .<= 1);

log10z = deepcopy(z);
log10z[pzind] .= log10.(z[pzind]);
log10z[nzind] .= -log10.(-z[nzind]);
log10z[zzind] .= 0;

fig = Figure(size = pres, fontsize = 54)
ax = Axis(
    fig[1, 1];
    title = "SMAB MARACOOS Glider Deployments (" * string(year(unix2datetime(NaNMath.minimum(gliderCTDarray[1].t)))) * "--" * string(year(unix2datetime(NaNMath.maximum(gliderCTDarray[end].t)))) * ")",
    xlabel = "Longitude",
    ylabel = "Latitude",
)
Makie.contourf!(x, y, log10z, colormap = :bukavu, levels = range(-6., 6., length = 128))
xlims!(lonmin, lonmax);
ylims!(latmin, latmax);

for ii = 1:length(gliderCTDarray)
    display(ii)
    Makie.scatter!(
        gliderCTDarray[ii].lon, 
        gliderCTDarray[ii].lat, 
        color = gliderCTDarray[ii].t, 
        colorrange = (minimum(gliderCTDarray[1].t), maximum(gliderCTDarray[end].t)), 
        colormap=:heat, 
        markersize=8
        )
end
fig
save("/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/SMAB_MARACOOS_glider.png", fig)
GLMakie.closeall()

#=pst = pst_electa;
glider1 = electaCTDraw;
glider2 = glider1;
#plot_glider_map(electaCTDraw, sylviaCTDraw, ps, pst);
include("slocumPlotMapTest.jl");
=#