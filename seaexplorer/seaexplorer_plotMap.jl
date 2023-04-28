# gong@vims.edu 2023-04-28
# this script plots glider data from SEA064's NORSE project

using GLMakie, NCDatasets, NaNMath

if (@isdefined jm) == false
    display("Loading NORSE SEA064 data.")
    include("run_seaexplorer.jl")
end

display("NORSE data loaded, begin plotting...")

pzlist = [0, -30, -60, -100, -150, -200, -300, -400, -600, -800];
pvarlist = ["ctemp", "saltA", "sigma0", "sndspd", "spice0"];

pvarlist = ["sndspd"];

for pvar in pvarlist
    for pz in pzlist
        if pvar == "ctemp"
            c1 = jm.ctemp;
            c2 = lbe.ctemp;
            cmin, cmax = -0.5, 9;
        elseif pvar == "saltA"
            c1 = jm.saltA;
            c2 = lbe.saltA;
            cmin, cmax = 34.5, 35.5;
        elseif pvar == "sigma0"
            c1 = jm.sigma0;
            c2 = lbe.sigma0;
            cmin, cmax = 27.3, 28.15;
        elseif pvar == "sndspd"
            c1 = jm.sndspd;
            c2 = lbe.sndspd;
            cmin, cmax = 1450, 1490;
        elseif pvar == "spice0"
            c1 = jm.spice0;
            c2 = lbe.spice0;
            cmin, cmax = -0.15, 1.2;
        end

        plotflag = pvar * "-" * string(abs(pz)) * "m";

        ptitle = "NORSE SEA064 " * pvar * " (" * string(abs(pz)) * "m)";
        pfname = "NORSE_SEA064_" * pvar * "_" * string(abs(pz)) * "m.png";
        zlo, zhi = pz-10, pz+10;

        x1 = jm.lon;
        y1 = jm.lat;
        z1 = jm.z;

        z1ind = findall(zlo .<= z1 .<= zhi);
        pind1 = z1ind;

        x2 = lbe.lon;
        y2 = lbe.lat;
        z2 = lbe.z;

        z2ind = findall(zlo .<= z2 .<= zhi);
        pind2 = z2ind;

        figoutdir = "/Users/gong/GitHub/jlglider/seaexplorer/figures/"
        bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
        bathyds = Dataset(bathypath,"r");

        lon = bathyds["lon"][:];
        lat = bathyds["lat"][:];

        # define plot boundaries
        latmin, latmax = 68, 73;
        lonmin, lonmax = -13, 17;
        tempmin, tempmax = -1, 9; # temperature

        #cmin, cmax = tempmin, tempmax;
        ms = 15;

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
            title = ptitle,
            xlabel = "Longitude",
            ylabel = "Latitude",
        )
        Makie.contourf!(x, y, z, xlims = (lonmin, lonmax), ylims = (latmin, latmax), levels = 128, colormap = :bukavu, colorrange = (-4000, 4000))
        Makie.scatter!(x1[pind1], y1[pind1], color = c1[pind1], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
        Makie.scatter!(x2[pind2], y2[pind2], color = c2[pind2], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
        Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = false)
        fig
        save(figoutdir * pfname, fig)
    end #pzlist
end #pvarlist
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
display("Plotting done.")