using GLMakie, NCDatasets, NaNMath

figoutdir = "/Users/gong/GitHub/jlglider/seaexplorer/figures/"
bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
bathyds = Dataset(bathypath,"r");

lon = bathyds["lon"][:];
lat = bathyds["lat"][:];

# define plot boundaries
latmin, latmax = 68, 73;
lonmin, lonmax = -13, 17;
tempmin, tempmax = -1, 9; # temperature

cmin, cmax = tempmin, tempmax;
ms = 10;

# approximate x-axis scaling to make it look "normal"
dlat = latmax - latmin;
dlon = lonmax - lonmin;
lat0 = 70;
xfac = sind(90-lat0);
yres = 1200;
pres = (ceil(yres * (xfac/(dlat/dlon))), yres);

# extract indices from the bathymetric data file
latind = findall(latmin-0.1 .<= lat .<= latmax+0.1);
lonind = findall(lonmin-0.1 .<= lon .<= lonmax+0.1);

z = bathyds["z"][lonind, latind];
zmin, zmax = -4000, 4000;
z[1,1] = -4000;
z[end,end] = 4000;

x = lon[lonind];
y = lat[latind];

fig = Figure(resolution = pres, fontsize = 32)
ax = Axis(
    fig[1, 1];
    title = "NORSE SEA064 surface temperature",
    xlabel = "Longitude",
    ylabel = "Latitude",
)
Makie.contourf!(x, y, z, xlims = (lonmin, lonmax), ylims = (latmin, latmax), levels = 128, colormap = :bukavu, colorrange = (-4000, 4000))
Makie.scatter!(x1[pind1], y1[pind1], color = c1[pind1], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
Makie.scatter!(x2[pind2], y2[pind2], color = c2[pind2], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = false)
fig
save(figoutdir * "NORSE_SEA064_ctemp.png", fig)

#= 
fig = Figure()

ga = GeoAxis(
    fig[1, 1]; # any cell of the figure's layout
    source = "+proj=latlong +datum=WGS84", # the CRS in which you want to plot
    dest = "+proj=wintri",
    lonlims=(-10, 10), latlims = (65, 80),
    coastlines = true # plot coastlines from Natural Earth, as a reference.
)
scatter!(ga, -120:15:120, -60:7.5:60; color = -60:7.5:60, strokecolor = (:black, 0.2))
fig
 =#