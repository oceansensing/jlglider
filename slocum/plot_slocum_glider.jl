using NaNMath, GibbsSeaWater, Dates, Interpolations
using GLMakie, ColorSchemes

import slocumFunc: datetick

# plot CTD data using Makie
# the plotting code will be refactored into a function of its own in the next revision


#function plot_glider_ctd(glider, figoutdir, ps)
    #if (@isdefined figoutdir) == false
        #figoutdir = "/Users/gong/Research/electa-20221103-passengers/figures/";
    #    rootdir = "/Users/gong/oceansensing Dropbox/C2PO/MARACOOS";
    #    figoutdir = rootdir * "/electa-20230320-maracoos/figures/";
    #    rootdir = "~/GitHub/jlglider/slocum/";
    #    figoutdir = figoutdir;
    #end
    mission = glider.mission;
    glider = glider.glider;

    #if (@isdefined ps) == false
    pint = ps.pint; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
    iday = ps.iday; # day intervals for plotting
    ms = ps.ms;
    tsms = ps.tsms;
    pres = ps.pres;
    tspres = ps.tspres;
    #end

    # setting x and y axes for plotting
    #td = dtctdf[1:pint:end];
    #xf = tctdf[1:pint:end]; 
    #yf = zzf[1:pint:end];

    #x = tctd;
    #xdt = unix2datetime.(tctd); 
    #yday = Dates.dayofyear.(xdt);
    #uyday = unique(yday);
    #hour = Dates.Hour.(xdt);
    #minute = Dates.Minute.(xdt);

    #tickind = Vector{Int64}(undef, length(uyday));
    #tickind = [];
    #for i = 1:length(uyday)
    #    t0 = findall((yday .== uyday[i]) .& (hour .== Hour(0)));
    #    if isempty(t0) == false
    #        push!(tickind, t0[1]);
    #    end
    #end
    #x = tctd;
    x = glider.t;
    xdt, xtick, xticklabel = datetick(x);
    #y = zzraw;
    y = glider.z;

    # plotting conservative temperature
    #z = ctempraw[1:pint:end];
    z = glider.ctemp;
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " Conservative Temperature",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (x[tickind], string.(Dates.Date.(xdt[tickind])));
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_ctemp.png", fig)

    #=
    # plotting absolute salinity
    z = saltAraw[1:pint:end];
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " Absolute Salinity",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xf[1]:86400*iday:xf[end], string.(Date.(td[1]:Day(iday):td[end])))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_saltA.png", fig)

    # plotting sigma0
    z = sigma0raw[1:pint:end];
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " Potential Density",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_sigma0.png", fig)

    # plotting spice0
    z = spice0raw[1:pint:end];
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " Spiciness0",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:balance, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_spice0.png", fig)

    # plotting sound speed
    z = sndspdraw[1:pint:end];
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " Sound Speed",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_soundspeed.png", fig)

    # plotting Chl-a
    x = chlatime;
    xdt, xtick, xticklabel = datetick(x);
    #xdt = chladtime;
    y = chlaz;
    z = chlaraw[1:pint:end,2];
    zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmax = 0.8;
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " Chl-a Concentration",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_chla.png", fig)

    # plotting CDOM
    x = cdomtime;
    xdt, xtick, xticklabel = datetick(x);
    y = cdomz;
    z = cdomraw[1:pint:end,2];
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = 0.0;
    zmax = 1.0;
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " CDOM",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_cdom.png", fig)

    # plotting BB700
    x = bb700time;
    xdt, xtick, xticklabel = datetick(x);
    y = bb700z;
    z = bb700raw[1:pint:end,2];
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = 0.0;
    zmax = 0.0004;
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " BB700",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glider * "_bb700.png", fig)

    # plotting BSI PAR
    x = bpartime;
    xdt, xtick, xticklabel = datetick(x);
    #xdt = bpardtime;
    y = bparz;
    z = log10.(bparraw[1:pint:end,2]);
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    zmin = -1.0;
    zmax = 4.0; 
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " PAR",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
    #ax.xticks = (x[1]:86400:x[end], string.(Date.(xdt[1]:Day(1):xdt[end])))
    ax.xticks = (xtick, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Log(PAR)")
    fig
    save(figoutdir * mission * "_" * glider * "_bsipar.png", fig)

    # plotting T/S diagram
    x = saltAraw[1:pint:end]; 
    y = ctempraw[1:pint:end];
    z = sigma0raw[1:pint:end];
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = tspres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * " T/S",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Sigma0")
    fig
    save(figoutdir * mission * "_" * glider * "_TS.png", fig)

    # plotting Density-Spiciness diagram
    x = spice0raw[1:pint:end]; 
    y = sigma0raw[1:pint:end];
    z = sndspdraw[1:pint:end];
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = tspres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glider * "Sigma0-Spice0",
        xlabel = "Spiciness",
        ylabel = "Potential Density Anomaly"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Sound Speed")
    fig
    save(figoutdir * mission * "_" * glider * "_sigma0spice0.png", fig)
    =#
#end