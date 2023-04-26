

figoutdir = "/Users/gong/GitHub/jlglider/microrider/figures/";

xfast = mrr.t_fast[:];
yfast = mrr.P_fast[:];
z1fast = mrr.sh1[:];
z2fast = mrr.sh2[:];

x = xfast;
y = yfast;
z = z1fast;

zmin = -0.02;
zmax = 0.02; 

pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 3; # day intervals for plotting
ms = 6;
pres = (1200, 800)

fig = Figure(resolution = pres)
ax = Axis(fig[1, 1],
    title = project * ": " * mission * " " * glider * " - dat_" * string(profileid, pad = 4) * " Shear 1",
    xlabel = "Time",
    ylabel = "Pressure"
)
Makie.lines!(x, y, color=z, colormap=:jet, linewidth=ms, colorrange=(zmin, zmax))
#ax.xticks = (xf[1]:86400*iday:xf[end], string.(Date.(td[1]:Day(iday):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
fig
save(figoutdir * project * "_" * mission * "_" * glider * "_" * string(profileid, pad = 4) * "_sh1.png", fig)