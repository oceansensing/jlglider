# 2025-03-17: DG

using NCDatasets, NaNMath, GLMakie

i = 2; # Electa's second mission during NESMA 2024
g = gliderCTDarray[i];

bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
bathyds = Dataset(bathypath,"r");

lon = bathyds["lon"][:];
lat = bathyds["lat"][:];

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
latind = findall(latmin-0.1 .<= lat .<= latmax+0.1);
lonind = findall(lonmin-0.1 .<= lon .<= lonmax+0.1);

z = Float64.(bathyds["z"][lonind, latind]); # recasting as Float64 to fix a StackOverFlow error seen in GLMakie 0.6.0.
x = lon[lonind];
y = lat[latind];

pzind = findall(z .> 1);
nzind = findall(z .< -1);
zzind = findall(-1 .<= z .<= 1);

logzflag = 0;
if logzflag == 1
    log10z = deepcopy(z);
    log10z[pzind] .= log10.(z[pzind]);
    log10z[nzind] .= -log10.(-z[nzind]);
    log10z[zzind] .= 0;
    zrange = (-4, 4);
    zp = log10z;
else
    zp = z;
    zrange = (-6000, 6000);
end

#plot(g.lon[gind], g.lat[gind])

varname = "ctemp"
plottitle = uppercase(gliderCTDarray[i].project) * " " * gliderCTDarray[i].glidertype * " " * varname * " (" * string(year(unix2datetime(NaNMath.minimum(gliderCTDarray[i].t)))) * "-" * string(year(unix2datetime(NaNMath.maximum(gliderCTDarray[i].t)))) * ")";
plotname = uppercase(gliderCTDarray[i].project) * "_" * lowercase(gliderCTDarray[i].glidername) * "_" * varname * ".png" 

fig = Figure(size = pres, fontsize = 64)
ax = Axis(
    fig[1, 1];
    title = plottitle,
    xlabel = "Longitude",
    ylabel = "Latitude",
)
Makie.contourf!(x, y, zp, colormap = ColorSchemes.bukavu, levels = range(zrange[1], zrange[2], length = 128))
xlims!(lonmin, lonmax);
ylims!(latmin, latmax);

for ii = 1:length(gliderCTDarray)
    zind = findall(minimum(pzrange) .<= gliderCTDarray[ii].z .<= maximum(pzrange));
    local x = gliderCTDarray[ii].lon[zind];
    local y = gliderCTDarray[ii].lat[zind];
    if varname == "saltA"
        local c = gliderCTDarray[ii].saltA[zind];
        cmin = pst[i].saltmin;
        cmax = pst[i].saltmax;
        cmap = ColorSchemes.haline;
    elseif varname == "ctemp"
        local c = gliderCTDarray[ii].ctemp[zind];
        cmin = pst[i].tempmin;
        cmax = pst[i].tempmax;
        cmap = ColorSchemes.thermal;
    elseif varname == "spice0"
        local c = gliderCTDarray[ii].spice0[zind];
        cmin = pst[i].spice0min;
        cmax = pst[i].spice0max;
        cmap = ColorSchemes.jet;
    end
    display(ii)
    GLMakie.scatter!(x, y, color=c, colorrange=(cmin,cmax), colormap=cmap, markersize=20);
    #Colorbar(fig[1, 2], limits = (cmin,cmax), colormap=cmap, flipaxis=true, label=varname)
end
fig
save(pst[i].figoutdir * plotname, fig)
GLMakie.closeall()