using Glider
include("/Users/gong/GitHub/jlglider/common/C2PO.jl")
using .C2PO: yearday2datetime, datetime2yearday
using NaNMath, GibbsSeaWater, Dates, Interpolations, Statistics, NCDatasets
using GLMakie, ColorSchemes

display("Plotting plotGliderTSarray...")
figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/NESMA-PASSENGERS/";

tspres = (1600, 1200);
fs = ps[1].fs;
ms = 1;
project = pst[end].project;
yyyy0 = string(Dates.year(unix2datetime(gliderCTDarray[1].t[1])), pad=4);
yyyyN = string(Dates.year(unix2datetime(gliderCTDarray[end].t[end])), pad=4);

saltmin = minimum([Float64(pst[i].saltmin) for i in 1:length(pst)]);
saltmax = maximum([Float64(pst[i].saltmax) for i in 1:length(pst)]);
tempmin = minimum([Float64(pst[i].tempmin) for i in 1:length(pst)]);
tempmax = maximum([Float64(pst[i].tempmax) for i in 1:length(pst)]);
spice0min = minimum([Float64(pst[i].spice0min) for i in 1:length(pst)]);
spice0max = maximum([Float64(pst[i].spice0max) for i in 1:length(pst)]);
sigma0min = minimum([Float64(pst[i].sigma0min) for i in 1:length(pst)]);
sigma0max = maximum([Float64(pst[i].sigma0max) for i in 1:length(pst)]);
sndspdmin = minimum([Float64(pst[i].sndspdmin) for i in 1:length(pst)]);
sndspdmax = maximum([Float64(pst[i].sndspdmax) for i in 1:length(pst)]); 
tmin = minimum(gliderCTDarray[1].t);
tmax = maximum(gliderCTDarray[end].t);

#c/spice0/sigma0 as color 3D
figSS = Figure(; size = tspres, fontsize = fs)
ax = Axis3(figSS[1, 1],
    title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " c/spice0/sigma0 - temp",
    xlabel = "Spice0",
    ylabel = "Sounds Speed (m/s)",
    zlabel = "Sigma0 (kg/m^3)"
)
xmin = spice0min;
xmax = spice0max;
ymin = spice0min;
ymax = spice0max;
zmin = pmin;
zmax = pmax;
cmin = tempmin;
cmax = tempmax;
xlims!(ax, xmin, xmax);
ylims!(ax, ymin, ymax);
zlims!(ax, zmin, zmax);

for i = 1:length(gliderCTDarray)
    gliderCTD = gliderCTDarray[i];
    x = gliderCTD.spice0;
    y = gliderCTD.sndspd;
    z = gliderCTD.sigma0;
    c = gliderCTD.ctemp;
    Makie.scatter!(x, y, z, color=c, colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
end
Colorbar(figSS[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="Conservative Temperature")
figSS
#save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_SS-Spice0-sigma0-ctemp.png", figSS);
#GLMakie.closeall()
