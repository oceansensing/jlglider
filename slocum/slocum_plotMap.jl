# gong@vims.edu 2023-04-28
# this script plots glider data from SEA064's NORSE project
workdir = "/Users/gong/GitHub/jlglider/slocum/"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

using GLMakie, NCDatasets, NaNMath, Dates, Interpolations
import seaexplorer_functions: seaexplorer_MR_laur_load

if (@isdefined jm) == false
    display("Loading PASSENGERS electa data.")
    include("run_passengers_glider_2023.jl")
end
display("PASSENGERS electa data loaded, begin plotting...")

pzlist = [0, -10, -20, -30, -40, -50, -60, -80, -100, -150, -200, -300, -400, -600, -800];
pvarlist = ["ctemp", "saltA", "sigma0", "sndspd", "spice0"];

#pzlist = [0];
pzlist = [0, -10, -20, -30, -40, -50, -60, -80, -90, -100, -125, -150, -175, -200, -225, -250, -275, -300, -325, -350, -375, -400, -450, -500, -550, -600, -650, -700, -750, -800, -850, -900];
pvarlist = ["spice0"]

region = "GS" # GS, MAB, or ALL

# define plot boundaries
if region == "GS"
    latmin, latmax = 37, 40;
    lonmin, lonmax = -65, -61;
elseif region == "MAB"
    latmin, latmax = 69, 71;
    lonmin, lonmax = 0, 10;
else
    latmin, latmax = 68, 73;
    lonmin, lonmax = -13, 17;
end

ms = 15;

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
            cmin, cmax = 27.3, 28.2;
        elseif pvar == "sndspd"
            c1 = jm.sndspd;
            c2 = lbe.sndspd;
            cmin, cmax = 1450, 1490;
        elseif pvar == "spice0"
            c1 = jm.spice0;
            c2 = lbe.spice0;
            cmin, cmax = -0.15, 1.2;
        elseif pvar == "epsilon"
            if (@isdefined eps1) == false
                #include("seaexplorer_plotMR.jl")
                lon1, lat1, lon2, lat2, eps1, eps2 = seaexplorer_MR_laur_load(jmpld1d, lbepld1d, pz, 10.0);
            end
            c1 = eps1;
            c2 = eps2;
            cmin, cmax = -11.5, -5.75;
        end

        plotflag = pvar * "-" * string(abs(pz)) * "m";

        ptitle = "NORSE SEA064 " * region * " " * pvar * " (" * string(abs(pz)) * "m)";
        
        if (region == "JM") | (region == "LBE")
            pfname = "NORSE_SEA064_" * region * "_" * pvar * "_" * string(abs(pz)) * "m.png";
        else
            pfname = "NORSE_SEA064_" * pvar * "_" * string(abs(pz)) * "m.png";
        end

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
        if pvar != "epsilon"
            if region != "LBE"
                Makie.scatter!(x1[pind1], y1[pind1], color = c1[pind1], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
            end
            if region != "JM"
                Makie.scatter!(x2[pind2], y2[pind2], color = c2[pind2], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
            end
        elseif pvar == "epsilon"
            c1 = []; c2 = [];
            log10eps1 = log10.(eps1);
            log10eps2 = log10.(eps2);
            gind1 = findall(isnan.(log10eps1) .!= true);
            gind2 = findall(isnan.(log10eps2) .!= true);
            c1 = log10eps1[gind1];
            c2 = log10eps2[gind2];
            if region != "LBE"
                Makie.scatter!(lon1[end,gind1], lat1[end,gind1], color = c1, colormap=:jet, markersize=ceil(ms*1.2), colorrange=(cmin, cmax), nan_color = RGBAf(0,0,0,0));
            end
            if region != "JM"
                Makie.scatter!(lon2[end,gind2], lat2[end,gind2], color = c2, colormap=:jet, markersize=ceil(ms*1.2), colorrange=(cmin, cmax), nan_color = RGBAf(0,0,0,0));
            end
        end
        Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = false)
        fig
        save(figoutdir * pfname, fig)
    end #pzlist
end #pvarlist

lon1, lat1, lon2, lat2, eps1, eps2 = seaexplorer_MR_laur_load(jmpld1d, lbepld1d, -10, 10.0);

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