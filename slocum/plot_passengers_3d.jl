# 2025-03-17: DG
#
# make 3D plot of NESMA-PASSENGERS electa deployment

using NCDatasets, NaNMath, GLMakie, ColorSchemes

i = 2; # Electa's second mission during NESMA 2024
g = gliderCTDarray[i];

bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
bathyds = Dataset(bathypath,"r");

lonb = bathyds["lon"][:];
latb = bathyds["lat"][:];

#latmin, latmax = 37, 40;
#lonmin, lonmax = -65.2, -59.8;
lonmin = pst[i].lonmin 
lonmax = pst[i].lonmax
latmin = pst[i].latmin
latmax = pst[i].latmax

# approximate x-axis scaling to make it look "normal"
dlat = latmax - latmin;
dlon = lonmax - lonmin;
lat0 = NaNMath.mean(gliderCTDarray[i].lat);
xfac = sind(90-lat0);
yres = 2000;
pres = (abs(ceil(yres * (xfac/(dlat/dlon)))), abs(yres));

# extract indices from the bathymetric data file
latbind = findall(latmin-0.1 .<= latb .<= latmax+0.1);
lonbind = findall(lonmin-0.1 .<= lonb .<= lonmax+0.1);

zb = Float64.(bathyds["z"][lonbind, latbind]); # recasting as Float64 to fix a StackOverFlow error seen in GLMakie 0.6.0.
xb = lonb[lonbind];
yb = latb[latbind];

pzbind = findall(zb .> 1);
nzbind = findall(zb .< -1);
zzbind = findall(-1 .<= zb .<= 1);

logzflag = 0;
if logzflag == 1
    log10zb = deepcopy(zb);
    log10zb[pzind] .= log10.(zb[pzbind]);
    log10zb[nzind] .= -log10.(-zb[nzbind]);
    log10zb[zzind] .= 0;
    zrange = (-4, 4);
    zbp = log10zb;
else
    zbp = zb;
    zrange = (-6000, 6000);
end

#plot(g.lon[gind], g.lat[gind])

varname = "ctemp"
plottitle = uppercase(gliderCTDarray[i].project) * " " * gliderCTDarray[i].glidername * " " * varname * " (" * Dates.format(unix2datetime(NaNMath.minimum(gliderCTDarray[i].t)), "yyyymmdd") * ")";
plotname = uppercase(gliderCTDarray[i].project) * "_" * lowercase(gliderCTDarray[i].glidername) * "_" * varname * "_" * Dates.format(unix2datetime(NaNMath.minimum(gliderCTDarray[i].t)), "yyyymmdd") * ".png" 

fig = Figure(size = pres, fontsize = 64)
ax = Axis3(
    fig[1, 1];
    title = plottitle,
    titlegap=0,
    xlabel = "",
    ylabel = "",
    zlabel = "",
    azimuth=π/6,    # 30° horizontal rotation
    elevation=π/6,   # 30° above xy-plane
    aspect=(1, 1, 0.3)
)

# Adjust layout to give more space at the top
rowsize!(fig.layout, 1, Relative(0.9))  # Reduce plot height to 90% of figure, leaving space for title

#GLMakie.contourf!(xb, yb, zbp, colormap = ColorSchemes.bukavu, levels = range(zrange[1], zrange[2], length = 128))
GLMakie.contourf!(xb, yb, zbp, colormap = ColorSchemes.bukavu, levels = range(zrange[1], zrange[2], length = 128), transformation=(:xy, -1000))
#GLMakie.surface!(xb, yb, zbp, colormap = ColorSchemes.topo)
xlims!(lonmin, lonmax);
ylims!(latmin, latmax);

pzrange = [-1000, 0];
zind = findall(minimum(pzrange) .<= gliderCTDarray[i].z .<= maximum(pzrange));
x = gliderCTDarray[i].lon[zind];
y = gliderCTDarray[i].lat[zind];
z = gliderCTDarray[i].z[zind];

if varname == "saltA"
    c = gliderCTDarray[i].saltA[zind];
    cmin = pst[i].saltmin;
    cmax = pst[i].saltmax;
    cmap = ColorSchemes.haline;
elseif varname == "ctemp"
    c = gliderCTDarray[i].ctemp[zind];
    cmin = pst[i].tempmin;
    cmax = pst[i].tempmax;
    cmap = ColorSchemes.thermal;
elseif varname == "spice0"
    c = gliderCTDarray[i].spice0[zind];
    cmin = pst[i].spice0min;
    cmax = pst[i].spice0max;
    cmap = ColorSchemes.jet;
end
display(i)
GLMakie.scatter!(x, y, z, color=c, colorrange=(cmin,cmax), colormap=ColorSchemes.jet, markersize=30);
#Colorbar(fig[1, 2], limits = (cmin,cmax), colormap=cmap, flipaxis=true, label=varname)

fig
save(pst[i].figoutdir * plotname, fig)
GLMakie.closeall()