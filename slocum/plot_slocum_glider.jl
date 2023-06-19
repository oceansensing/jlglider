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
    mission = gliderCTD.mission;
    glidername = gliderCTD.glidername;

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
    x0 = datetime2unix.(DateTime("2022-01-01"));
    x = gliderCTD.t;
    xdt, xtick, xticklabel = datetick(x);
    #y = zzraw;
    y = gliderCTD.z;

    # plotting conservative temperature
    z = gliderCTD.ctemp;
    zmin = NaNMath.minimum(z);
    #zmin = 18;
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glidername * " Conservative Temperature",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick .- x0, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glidername * "_ctemp.png", fig)

    # plotting absolute salinity
    #z = saltAraw[1:pint:end];
    z = gliderCTD.saltA;
    zmin = NaNMath.minimum(z);
    #zmin = 36;
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glidername * " Absolute Salinity",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick .- x0, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glidername * "_saltA.png", fig)

    # plotting sigma0
    #z = sigma0raw[1:pint:end];
    z = gliderCTD.sigma0;
    zmin = NaNMath.minimum(z);
    #zmin = 24.0;
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glidername * " Potential Density",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick .- x0, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glidername * "_sigma0.png", fig)

    # plotting spice0
    #z = spice0raw[1:pint:end];
    z = gliderCTD.spice0;
    zmin = NaNMath.minimum(z);
    #zmin = 4.0;
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glidername * " Spiciness0",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x .- x0, y, color=z, colormap=:balance, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick .- x0, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :balance, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glidername * "_spice0.png", fig)

    # plotting sound speed
    #z = sndspdraw[1:pint:end];
    z = gliderCTD.sndspd;
    #zmin = NaNMath.minimum(z);
    zmin = 1520;
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = pres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glidername * " Sound Speed",
        xlabel = "Time",
        ylabel = "Depth"
    )
    Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    ax.xticks = (xtick .- x0, xticklabel);
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
    fig
    save(figoutdir * mission * "_" * glidername * "_soundspeed.png", fig)

    # plotting T/S diagram
    #x = saltAraw[1:pint:end]; 
    #y = ctempraw[1:pint:end];
    #z = sigma0raw[1:pint:end];
    x = gliderCTD.saltA;
    y = gliderCTD.ctemp;
    z = gliderCTD.sigma0;
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = tspres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glidername * " T/S",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Sigma0")
    fig
    save(figoutdir * mission * "_" * glidername * "_TS.png", fig)

    # plotting Density-Spiciness diagram
    #x = spice0raw[1:pint:end]; 
    #y = sigma0raw[1:pint:end];
    #z = sndspdraw[1:pint:end];
    x = gliderCTD.spice0;
    y = gliderCTD.sigma0;
    z = gliderCTD.sndspd;
    zmin = NaNMath.minimum(z);
    zmax = NaNMath.maximum(z); 
    fig = Figure(resolution = tspres)
    ax = Axis(fig[1, 1],
        title = mission * " " * glidername * "Sigma0-Spice0",
        xlabel = "Spiciness",
        ylabel = "Potential Density Anomaly"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Sound Speed")
    fig
    save(figoutdir * mission * "_" * glidername * "_sigma0spice0.png", fig)

#=    
    # plotting Chl-a
    if isdefined(gliderCHLA, :t) == true
        x0 = datetime2unix.(DateTime("2022-01-01"));
        #x = chlatime;
        x = gliderCHLA.t;
        xdt, xtick, xticklabel = datetick(x);
        y = gliderCHLA.z;
        z = gliderCHLA.var;
        zmin = NaNMath.minimum(z);
        #zmax = NaNMath.maximum(z); 
        zmax = 0.8;
        fig = Figure(resolution = (1200, 800))
        ax = Axis(fig[1, 1],
            title = gliderCHLA.mission * " " * gliderCHLA.glidername * " Chl-a Concentration",
            xlabel = "Time",
            ylabel = "Depth"
        )
        Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
        ax.xticks = (xtick .- x0, xticklabel);
        Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
        fig
        save(figoutdir * gliderCHLA.mission * "_" * gliderCHLA.glidername * "_chla.png", fig);
    end

    # plotting CDOM
    if isdefined(gliderCDOM, :t) == true
        x0 = datetime2unix.(DateTime("2022-01-01"));
        x = gliderCDOM.t;
        xdt, xtick, xticklabel = datetick(x);
        y = gliderCDOM.z;
        z = gliderCDOM.var;
        #zmin = NaNMath.minimum(z);
        #zmax = NaNMath.maximum(z); 
        zmin = 0.0;
        zmax = 1.0;
        fig = Figure(resolution = (1200, 800))
        ax = Axis(fig[1, 1],
            title = gliderCDOM.mission * " " * gliderCDOM.glidername * " CDOM",
            xlabel = "Time",
            ylabel = "Depth"
        )
        Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
        ax.xticks = (xtick .- x0, xticklabel);
        Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
        fig
        save(figoutdir * gliderCDOM.mission * "_" * gliderCDOM.glidername * "_cdom.png", fig);
    end

    # plotting BB700
    if isdefined(gliderBB700, :t) == true
        x0 = datetime2unix.(DateTime("2022-01-01"));
        x = gliderBB700.t;
        xdt, xtick, xticklabel = datetick(x);
        y = gliderBB700.z;
        z = gliderBB700.var;
        #zmin = NaNMath.minimum(z);
        #zmax = NaNMath.maximum(z); 
        zmin = 0.0;
        zmax = 0.0004;
        fig = Figure(resolution = (1200, 800))
        ax = Axis(fig[1, 1],
            title = gliderBB700.mission * " " * gliderBB700.glidername * " BB700",
            xlabel = "Time",
            ylabel = "Depth"
        )
        Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
        ax.xticks = (xtick .- x0, xticklabel);
        Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
        fig
        save(figoutdir * gliderBB700.mission * "_" * gliderBB700.glidername * "_bb700.png", fig);
    end

    # plotting BSI PAR
    if isdefined(gliderBPAR, :t) == true
        x0 = datetime2unix.(DateTime("2022-01-01"));
        x = gliderBPAR.t;
        xdt, xtick, xticklabel = datetick(x);
        y = gliderBPAR.z;
        z = log10.(gliderBPAR.var);
        #zmin = NaNMath.minimum(z);
        #zmax = NaNMath.maximum(z);
        zmin = -1.0;
        zmax = 4.0; 
        fig = Figure(resolution = (1200, 800))
        ax = Axis(fig[1, 1],
            title = gliderBPAR.mission * " " * gliderBPAR.glidername * " PAR",
            xlabel = "Time",
            ylabel = "Depth"
        )
        Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=6, colorrange=(zmin, zmax))
        #ax.xticks = (x[1]:86400:x[end], string.(Date.(xdt[1]:Day(1):xdt[end])))
        ax.xticks = (xtick .- x0, xticklabel);
        Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false, label="Log(PAR)")
        fig
        save(figoutdir * gliderBPAR.mission * "_" * gliderBPAR.glidername * "_bsipar.png", fig);
    end

#end
=#