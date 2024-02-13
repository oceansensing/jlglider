using CSV, DataFrames
using NaNMath, GibbsSeaWater, Dates, Interpolations
using GLMakie, ColorSchemes

include("seaexplorerType.jl")
include("seaexplorerFunc.jl")
import .seaexplorerType: plotSetting, plotStruct, engStruct, ctdStruct, sciStruct
import .seaexplorerFunc: datetick, yearday

glidername = "SEA064"
projectname = "NORSE"
missionid = ["37"; "38"; "48"];
missiondate = ["20221021"; "20221102"; "20231112"];
region = ["janmayen"; "lofoten"; "janmayen"];

for i = 1:3
    if i == 1
        gdata = jm22;
    elseif i == 2
        gdata = lbe22;
    elseif i == 3
        gdata = jm23;
    end

    gliderDataDir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
    gliderMissionDir = lowercase(glidername) * "-" * missiondate[i] * "-" * lowercase(projectname) * "-" * region[i] * "-complete";
    gliderADCPdir = gliderDataDir * gliderMissionDir * "/ad2cp/" * "m" * missionid[i] * "_processed/";
    gliderADCPpath = gliderADCPdir * "abs_dive_ocean_vel.csv";
    figoutdir = gliderADCPdir * "figures/";

    #        df = CSV.read(navfilepath, header=1, delim=";", DataFrame, buffer_in_memory=true);

    gdf = CSV.read(gliderADCPpath, header=1, delim=",", DataFrame);

    yo = Int64.(gdf.yo_number);
    t = collect(gdf.time_midpoint);
    u = collect(gdf.u_ocean_vel);
    v = collect(gdf.v_ocean_vel);
    z = collect(gdf.depth_bins);

    #t0 = datetime2unix.(DateTime("2022-01-01"));
    xdt, xtick, xticklabel = datetick(t);

    tdt = unix2datetime.(t); 
    yday = yearday.(tdt);

    pres = (1200, 800);
    fs = 24;
    ms = 12;
    zmin, zmax = -0.5, 0.5;

    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/Presentations/OSM2024/";

    figU = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figU[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " E-W ADCP velocity (post-processed)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yday, z, color=u, colormap=:balance, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figU[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = true, label = "u (m/s)")
    fig
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_U.png", figU)
    GLMakie.closeall()

    figV = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figV[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " N-S ADCP velocity (post-processed)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yday, z, color=v, colormap=:balance, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figV[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = true, label = "v (m/s)")
    fig
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_V.png", figV)
    GLMakie.closeall()

    figUrt = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figUrt[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " E-W ADCP velocity (realtime)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yearday.(gdata.t), gdata.z, color=gdata.ad2cp_Ueast, colormap=:balance, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figUrt[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = true, label = "u (m/s)")
    fig
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_Urt.png", figUrt)
    GLMakie.closeall()

    figVrt = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figVrt[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " N-S ADCP velocity (realtime)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yearday.(gdata.t), gdata.z, color=gdata.ad2cp_Unorth, colormap=:balance, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figVrt[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = true, label = "v (m/s)")
    fig
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_Vrt.png", figVrt)
    GLMakie.closeall()
end