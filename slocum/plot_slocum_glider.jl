using NaNMath, GibbsSeaWater, Dates, Interpolations
using GLMakie, ColorSchemes

# plot CTD data using Makie
# the plotting code will be refactored into a function of its own in the next revision

if (@isdefined figoutdir) == false
    #figoutdir = "/Users/gong/Research/electa-20221103-passengers/figures/";
    rootdir = "/Users/gong/oceansensing Dropbox/C2PO/MARACOOS";
    figoutdir = rootdir * "/electa-20230320-maracoos/figures/";
    mission = "MARACOOS";
    glider = "electa";
end

pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 3; # day intervals for plotting

# setting x and y axes for plotting
td = dtctdf[1:pint:end];
xf = tctdf[1:pint:end]; 
yf = zzf[1:pint:end];

x = tctd;
y = zzraw;

# plotting conservative temperature
z = ctempraw[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = mission * " " * glider * " Conservative Temperature",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(x, y, color=z, colormap=:thermal, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (xf[1]:86400*iday:xf[end], string.(DateTime.(td[1]:Day(iday):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :thermal, flipaxis = false)
fig
save(figoutdir * mission * "_" * glider * "_ctemp.png", fig)

# plotting absolute salinity
z = saltAraw[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = mission * " " * glider * " Absolute Salinity",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(x, y, color=z, colormap=:haline, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (xf[1]:86400*iday:xf[end], string.(Date.(td[1]:Day(iday):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :haline, flipaxis = false)
fig
save(figoutdir * mission * "_" * glider * "_saltA.png", fig)

# plotting sigma0
z = sigma0raw[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = mission * " " * glider * " Potential Density",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(x, y, color=z, colormap=:dense, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (xf[1]:86400*iday:xf[end], string.(Date.(td[1]:Day(iday):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :dense, flipaxis = false)
fig
save(figoutdir * mission * "_" * glider * "_sigma0.png", fig)

# plotting spice0
z = spice0raw[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = mission * " " * glider * " Spiciness0",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(x, y, color=z, colormap=:balance, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (xf[1]:86400*iday:xf[end], string.(Date.(td[1]:Day(iday):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = false)
fig
save(figoutdir * mission * "_" * glider * "_spice0.png", fig)

# plotting sound speed
z = sndspdraw[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = mission * " " * glider * " Sound Speed",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(x, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (xf[1]:86400*iday:xf[end], string.(Date.(td[1]:Day(iday):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
fig
save(figoutdir * mission * "_" * glider * "_soundspeed.png", fig)

## plotting Chl-a
#z = chlaraw[1:pint:end,2];
#zmin = NaNMath.minimum(z);
#zmax = NaNMath.maximum(z); 
#fig = Figure(resolution = (1200, 800))
#ax = Axis(fig[1, 1],
#    title = mission * " " * glider * " Chl-a Concentration",
#    xlabel = "Time",
#    ylabel = "Depth"
#)
#Makie.scatter!(x, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
#ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
#Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
#fig
#save(figoutdir * mission * "_" * glider * "_chla.png", fig)

# plotting T/S diagram
x = saltAraw[1:pint:end]; 
y = ctempraw[1:pint:end];
z = sigma0raw[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1000, 1000))
ax = Axis(fig[1, 1],
    title = mission * " " * glider * " T/S",
    xlabel = "Absolute Salinity",
    ylabel = "Conservative Temperature"
)
Makie.scatter!(x, y, color=z, colormap=:jet, markersize=2, colorrange=(zmin, zmax))
#ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Sigma0")
fig
save(figoutdir * mission * "_" * glider * "_TS.png", fig)
