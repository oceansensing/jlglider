using Glider
include("/Users/gong/GitHub/jlglider/common/C2PO.jl")
using .C2PO: yearday2datetime, datetime2yearday
using NaNMath, GibbsSeaWater, Dates, Interpolations, Statistics, NCDatasets
using GLMakie, ColorSchemes

zflag = "snds"

display("Plotting plotGliderTSarray...")
figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/NESMA-PASSENGERS/";

tspres = (1600, 1200);
fs = ps[1].fs;
ms = 1;
project = pst[end].project;
yyyy0 = string(Dates.year(unix2datetime(gliderCTDarray[1].t[1])), pad=4);
yyyyN = string(Dates.year(unix2datetime(gliderCTDarray[end].t[end])), pad=4);

saltmin = minimum([Float64(pst[i].saltmin) for i in 1:length(pst)]);
saltmax = maximum([Float64(pst[i].saltmax) for i in 1:length(pst)]);
tempmin = minimum([Float64(pst[i].tempmin) for i in 1:length(pst)]);
tempmax = maximum([Float64(pst[i].tempmax) for i in 1:length(pst)]);
spice0min = minimum([Float64(pst[i].spice0min) for i in 1:length(pst)]);
spice0max = maximum([Float64(pst[i].spice0max) for i in 1:length(pst)]);
sigma0min = minimum([Float64(pst[i].sigma0min) for i in 1:length(pst)]);
sigma0max = maximum([Float64(pst[i].sigma0max) for i in 1:length(pst)]);
sndspdmin = minimum([Float64(pst[i].sndspdmin) for i in 1:length(pst)]);
sndspdmax = maximum([Float64(pst[i].sndspdmax) for i in 1:length(pst)]); 
pmin = minimum([Float64(pst[i].pmin) for i in 1:length(pst)]);
pmax = maximum([Float64(pst[i].pmax) for i in 1:length(pst)]);
tmin = minimum(gliderCTDarray[1].t);
tmax = maximum(gliderCTDarray[end].t);

if lowercase(zflag[1:4]) == "sigm"
    #ctemp/saltA/sigma0 - c color 3D
    figTS = Figure(; size = tspres, fontsize = fs)
    ax = Axis3(figTS[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " T/S/sigma0 - c",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature",
        zlabel = "Sigma0 (kg/m^3)"
    )
    xmin = saltmin;
    xmax = saltmax;
    ymin = tempmin;
    ymax = tempmax;
    zmin = sigma0min;
    zmax = sigma0max;
    cmin = sndspdmin;
    cmax = sndspdmax;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);
    zlims!(ax, zmin, zmax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.saltA;
        y = gliderCTD.ctemp;
        z = gliderCTD.sigma0;
        c = gliderCTD.sndspd;
        Makie.scatter!(x, y, z, color=c, colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
    end
    Colorbar(figTS[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="Sound Speed (m/s)")
    figTS
    #save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_SA-TS-sigma0-c.png", figTS);
    #GLMakie.closeall()

elseif lowercase(zflag[1:4]) == "snds" 

    #ctemp/saltA/sigma0 - c color 3D
    figTS2 = Figure(; size = tspres, fontsize = fs)
    ax = Axis3(figTS2[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " T/S/c - sigma0",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature",
        zlabel = "Sound Speed (m/s)"
    )
    xmin = saltmin;
    xmax = saltmax;
    ymin = tempmin;
    ymax = tempmax;
    cmin = sigma0min;
    cmax = sigma0max;
    zmin = sndspdmin;
    zmax = sndspdmax;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);
    zlims!(ax, zmin, zmax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.saltA;
        y = gliderCTD.ctemp;
        c = gliderCTD.sigma0;
        z = gliderCTD.sndspd;
        Makie.scatter!(x, y, z, color=c, colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
    end
    Colorbar(figTS2[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="Sigma0 (kg/m^3)")
    figTS2
    #save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_SS-Spice0-sigma0-ctemp.png", figTS2);
    #GLMakie.closeall()
elseif lowercase(zflag[1:4]) == "pres" 

    #ctemp/press/sndspd - saltA color 3D
    figTS3 = Figure(; size = tspres, fontsize = fs)
    ax = Axis3(figTS3[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " CT/SA/P - c",
        xlabel = "Absolute Salinity (ppt)",
        ylabel = "Conservative Temperature (C)",
        zlabel = "Pressure (dbar)"
    )
    xmin = saltmin;
    xmax = saltmax;
    ymin = tempmin;
    ymax = tempmax;
    cmin = sndspdmin;
    cmax = sndspdmax;
    zmin = pmin;
    zmax = pmax;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);
    zlims!(ax, zmin, zmax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.saltA;
        y = gliderCTD.ctemp;
        c = gliderCTD.sndspd;
        z = gliderCTD.p*10;
        Makie.scatter!(x, y, z, color=c, colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
    end
    Colorbar(figTS3[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="c (m/s)")
    figTS3
    #save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_SA-CT-P-c.png", figTS3);
    #GLMakie.closeall()
end