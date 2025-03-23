workdir = "/Users/gong/GitHub/jlglider/slocum/"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

using GLMakie, NCDatasets, NaNMath, Dates, Interpolations, ColorSchemes
#import seaexplorer_functions: seaexplorer_MR_laur_load

region = "GS"

#lonrange = [-9.5 -6.0];
#latrange = [70.5 71.5]; 

lonrange = [-65.5 -61.0];
latrange = [38.0 40.0]; 

latmin, latmax = latrange[1], latrange[2];
lonmin, lonmax = lonrange[1], lonrange[2];

pzlist = [0, -10, -20, -30, -40, -50, -60, -80, -100, -150, -200];
pvarlist = ["ctemp", "saltA", "sigma0", "sndspd", "spice0"];

#pz = -50;
#pvar = "sndspd"
for pvar in pvarlist
    for pz in pzlist
        if pvar == "ctemp"
            c1 = glider1.ctemp;
            c2 = glider2.ctemp;
            cmin, cmax = pst.tempmin, pst.tempmax;
        elseif pvar == "saltA"
            c1 = glider1.saltA;
            c2 = glider2.saltA;
            cmin, cmax = pst.saltmin, pst.saltmax;
        elseif pvar == "sigma0"
            c1 = glider1.sigma0;
            c2 = glider2.sigma0;
            cmin, cmax = pst.sigma0min, pst.sigma0max;
        elseif pvar == "sndspd"
            c1 = glider1.sndspd;
            c2 = glider2.sndspd;
            cmin, cmax = pst.sndspdmin, pst.sndspdmax;
        elseif pvar == "spice0"
            c1 = glider1.spice0;
            c2 = glider2.spice0;
            cmin, cmax = pst.spice0min, pst.spice0max;
        elseif pvar == "epsilon"
            if (@isdefined eps1) == false
                #include("seaexplorer_plotMR.jl")
            #    lon1, lat1, lon2, lat2, eps1, eps2 = seaexplorer_MR_laur_load(jmpld1d, lbepld1d, pz, 10.0);
            end
            c1 = eps1;
            c2 = eps2;
            cmin, cmax = -11.5, -5.75;
        end

        plotflag = pvar * "-" * string(abs(pz)) * "m";

        ptitle = pst.mission * " " * pst.glidername * " " * pvar * " (" * string(abs(pz)) * "m)";

        if region in ["JM", "LBE", "GS", "MAB"]
            pfname = pst.mission * "_" * region * "_" * pvar * "_" * string(abs(pz)) * "m.png";
        else
            pfname = pst.mission * "_" * pvar * "_" * string(abs(pz)) * "m.png";
        end

        x1 = glider1.lon;
        y1 = glider1.lat;
        z1 = glider1.z;

        x2 = glider2.lon;
        y2 = glider2.lat;
        z2 = glider2.z;

        zlo, zhi = pz-5, pz+5;
        z1ind = findall(zlo .<= z1 .<= zhi);
        pind1 = z1ind;
        z2ind = findall(zlo .<= z2 .<= zhi);
        pind2 = z2ind;

        figoutdir = pst.figoutdir;
        bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
        bathyds = Dataset(bathypath,"r");

        lon = bathyds["lon"][:];
        lat = bathyds["lat"][:];

        # approximate x-axis scaling to make it look "normal"
        dlat = latmax - latmin;
        dlon = lonmax - lonmin;
        lat0 = NaNMath.mean([y1; y2]);
        xfac = sind(90-lat0);
        yres = 1200;
        pres = (ceil(yres * (xfac/(dlat/dlon))), yres);

        latind = findall(latmin-0.1 .<= lat .<= latmax+0.1);
        lonind = findall(lonmin-0.1 .<= lon .<= lonmax+0.1);
        x = lon[lonind];
        y = lat[latind];
        z = Float64.(bathyds["z"][lonind, latind]);

        zmin, zmax = -5400, 5400;
        z[1,1] = -5400;
        z[end,end] = 5400;

        fig = Figure(resolution = pres, fontsize = 36)
        ax = Axis(
            fig[1, 1];
            title = ptitle,
            xlabel = "Longitude",
            ylabel = "Latitude",
        )
        ms = 12
        #GLMakie.contourf!(x, y, z, xlims = (lonmin, lonmax), ylims = (latmin, latmax), levels = 128, colormap = :bukavu, colorrange = (-5400, 5400))
        GLMakie.contourf!(x, y, z, levels = 128, colormap = :bukavu)
        GLMakie.xlims!(ax, lonmin, lonmax);
        GLMakie.ylims!(ax, latmin, latmax);
        GLMakie.scatter!(x1[pind1], y1[pind1], color = c1[pind1], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
        GLMakie.scatter!(x2[pind2], y2[pind2], color = c2[pind2], colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
        Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = false)
        fig
        save(figoutdir * pfname, fig)
    end
end