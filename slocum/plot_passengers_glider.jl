using NaNMath, GibbsSeaWater, Dates, Interpolations
using GLMakie, ColorSchemes

# plot CTD data
pint = 1;
#hctemp = scatter(dtctd[1:pint:end], z[1:pint:end], zcolor=ctemp[1:pint:end], legend=false, title="PASSENGERS Glider Electa Conservative Temperature", clims=(18, 26), size = (1200, 800), colorbar = true, markersize = 2.5, markerstrokewidth = -1, framestyle=:box, seriescolor=:thermal)
#hsaltA = scatter(dtctd[1:pint:end], z[1:pint:end], zcolor=saltA[1:pint:end], legend=false, title="PASSENGERS Glider Electa Absolute Salinity", clims=(36, 37.2), size = (1200, 800), colorbar = true, markersize = 2.5, markerstrokewidth = -1, framestyle=:box, seriescolor=:haline)
#Plots.savefig(hctemp, figoutdir * "PASSENGERS_electa_ctemp.png");
#Plots.savefig(hsaltA, figoutdir * "PASSENGERS_electa_saltA.png");

td = dtctd[1:pint:end];
t = tctd[1:pint:end]; 
y = zz[1:pint:end];

# plotting conservative temperature
z = ctemp[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = "PASSENGERS VIMS Glider Electa Conservative Temperature",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(t, y, color=z, colormap=:thermal, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :thermal, flipaxis = false)
fig
save(figoutdir * "PASSENGERS_electa_ctemp.png", fig)

# plotting absolute salinity
z = saltA[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = "PASSENGERS VIMS Glider Electa Absolute Salinity",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(t, y, color=z, colormap=:haline, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :haline, flipaxis = false)
fig
save(figoutdir * "PASSENGERS_electa_saltA.png", fig)

# plotting sigma0
z = sigma0[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = "PASSENGERS VIMS Glider Electa Potential Density",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(t, y, color=z, colormap=:dense, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :dense, flipaxis = false)
fig
save(figoutdir * "PASSENGERS_electa_sigma0.png", fig)

# plotting spice0
z = spice0[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = "PASSENGERS VIMS Glider Electa Spiciness0",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(t, y, color=z, colormap=:balance, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = false)
fig
save(figoutdir * "PASSENGERS_electa_spice0.png", fig)

# plotting sound speed
z = sndspd[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    title = "PASSENGERS VIMS Glider Electa Sound Speed",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(t, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
fig
save(figoutdir * "PASSENGERS_electa_soundspeed.png", fig)


# plotting sound speed
x = saltA[1:pint:end]; 
y = ctemp[1:pint:end];
z = sndspd[1:pint:end];
zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = (1000, 1000))
ax = Axis(fig[1, 1],
    title = "PASSENGERS VIMS Glider Electa T/S",
    xlabel = "Absolute Salinity",
    ylabel = "Conservative Temperature"
)
Makie.scatter!(x, y, color=z, colormap=:jet, markersize=2, colorrange=(zmin, zmax))
#ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Sound Speed")
fig
save(figoutdir * "PASSENGERS_electa_TS.png", fig)
