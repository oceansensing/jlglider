module gliderPlot

using Glider
#include("/Users/gong/GitHub/jlglider/slocum/gliderPlotType.jl")
include("/Users/gong/GitHub/ocean_julia/C2PO.jl")
using .C2PO: yearday2datetime, datetime2yearday
#using .gliderType: plotSetting, plotStruct, ctdStruct, sciStruct
using NaNMath, GibbsSeaWater, Dates, Interpolations, Statistics, NCDatasets
using GLMakie, ColorSchemes

function datetick(unix_t::Array{Float64})
    x = unix_t
    xdt = unix2datetime.(unix_t); 
    yday = Dates.dayofyear.(xdt);
    uyday = unique(yday);
    hour = Dates.Hour.(xdt);
    minute = Dates.Minute.(xdt);
    #df = DateFormat("y-m-d");

    #tickind = Vector{Int64}(undef, length(uyday));
    tickind = [];
    for i = 1:length(uyday)
        t0 = findall((yday .== uyday[i]) .& (hour .== Hour(0)));
        if isempty(t0) == false
            push!(tickind, t0[1]);
        end
    end
    xtick = x[tickind];
    #xticklabel = string.(Dates.Date.(xdt[tickind]));
    xticklabel = [x[6:10] for x in string.(xdt[tickind])];
    return xdt, xtick, xticklabel    
end

function plot_glider_ctd(gliderCTD, ps, pst)
    #if (@isdefined figoutdir) == false
        #figoutdir = "/Users/gong/Research/electa-20221103-passengers/figures/";
    #    rootdir = "/Users/gong/oceansensing Dropbox/C2PO/MARACOOS";
    #    figoutdir = rootdir * "/electa-20230320-maracoos/figures/";
    #    rootdir = "~/GitHub/jlglider/slocum/";
    #    figoutdir = figoutdir;
    #end
    project = pst.project;
    glidername = pst.glidername;

    #if (@isdefined ps) == false
    pint = ps.pint; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
    iday = ps.iday; # day intervals for plotting
    ms = ps.ms;
    tsms = ps.tsms;
    pres = ps.pres;
    tspres = ps.tspres;
    fs = ps.fs;
    figoutdir = pst.figoutdir;

    nsig = 2.5; # number of standard deviations for color range
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
    #x0 = datetime2unix.(DateTime(string(Dates.year(unix2datetime(gliderCTD.t[1]))) * "-01-01"));
    yyyy0 = string(Dates.year(unix2datetime(gliderCTD.t[1])), pad=4);
    mm0 = string(Dates.month(unix2datetime(gliderCTD.t[1])), pad=2);
    dd0 = string(Dates.day(unix2datetime(gliderCTD.t[1])), pad=2);
    tdt = unix2datetime.(gliderCTD.t); 
    yday = datetime2yearday.(tdt);

    x = yday;
    xdt, xtick, xticklabel = datetick(x);
    #y = zzraw;
    y = gliderCTD.z;

    # plotting temperature
    z = gliderCTD.temp;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    zmin = pst.tempmin;
    zmax = pst.tempmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Temperature",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #=
    ax.xticks = (xtick .- x0, xticklabel);
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#

    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Temperature (°C)")
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_temp.png", fig)

    # plotting conductivity
    z = gliderCTD.cond;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    zmin = pst.condmin;
    zmax = pst.condmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Conductivity",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Conductivity",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Conductivity (S/m)")
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_cond.png", fig)

    # plotting salinity
    z = gliderCTD.salt;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = pst.saltmin;
    zmax = pst.saltmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Salinity",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Salinity",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Salinity (psu)")
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_salt.png", fig)

    # plotting conservative temperature
    z = gliderCTD.ctemp;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = pst.tempmin;
    zmax = pst.tempmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Conservative Temperature",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Conservative Temperature",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Conservative Temperature (°C)")
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_ctemp.png", fig)

    # plotting absolute salinity
    #z = saltAraw[1:pint:end];
    z = gliderCTD.saltA;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = pst.saltmin;
    zmax = pst.saltmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Absolute Salinity",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Absolute Salinity",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Absolute Salinity (g/kg)")
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_saltA.png", fig)

    # plotting sigma0
    #z = sigma0raw[1:pint:end];
    z = gliderCTD.sigma0;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    zmin = pst.sigma0min;
    zmax = pst.sigma0max;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = project * " " * uppercasefirst(glidername) * " Potential Density",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Potential Density",

        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Potential Density (kg/m^3)")
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_sigma0.png", fig)

    # plotting spice0
    #z = spice0raw[1:pint:end];
    z = gliderCTD.spice0;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    zmin = pst.spice0min;
    zmax = pst.spice0max;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Spiciness0",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Spice0",

        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Spice0")
    fig
    #save(figoutdir * project * "_" * glidername * "_spice0.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_spice0.png", fig)

    # plotting sound speed
    #z = sndspdraw[1:pint:end];
    z = gliderCTD.sndspd;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = pst.sndspdmin;
    zmax = pst.sndspdmax;

    fig = Figure(; size = pres, fontsize = fs);
    ax = Axis(fig[1, 1],
        #title = mission * " – " * uppercasefirst(glidername) * " sound speed",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Sound Speed",
        #title =  uppercasefirst(glidername) * " Sound Speed",
        xlabel = "Year Day",
        ylabel = "Depth (m)",
        #label = :bold
    )
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Sound Speed (m/s)")
    fig
    #save(figoutdir * project * "_" * glidername * "_soundspeed.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_soundspeed.png", fig)


    # plotting CT/SA diagram
    #x = saltAraw[1:pint:end]; 
    #y = ctempraw[1:pint:end];
    #z = sigma0raw[1:pint:end];
    x = gliderCTD.saltA;
    y = gliderCTD.ctemp;
    z = gliderCTD.sigma0;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    #xmin = NaNMath.mean(x) .- (nsig+0.5) * NaNMath.std(x);
    #xmax = NaNMath.mean(x) .+ (nsig+0.5) * NaNMath.std(x);
    #ymin = NaNMath.mean(y) .- (nsig+0.5) * NaNMath.std(y);
    #ymax = NaNMath.mean(y) .+ (nsig+0.5) * NaNMath.std(y); 
    #zmin = NaNMath.mean(z) .- nsig * NaNMath.std(z);
    #zmax = NaNMath.mean(z) .+ nsig * NaNMath.std(z);

    xmin = pst.saltmin;
    xmax = pst.saltmax;
    ymin = pst.tempmin;
    ymax = pst.tempmax;
    zmin = pst.sigma0min;
    zmax = pst.sigma0max;

    fig = Figure(; size = tspres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " CT/SA",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " CT/SA",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature"
    )
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label="Sigma0")
    fig
    #save(figoutdir * project * "_" * glidername * "_CT-SA.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_CA-SA.png", fig)

    # plotting Density-Spiciness diagram
    #x = spice0raw[1:pint:end]; 
    #y = sigma0raw[1:pint:end];
    #z = sndspdraw[1:pint:end];
    x = gliderCTD.spice0;
    y = gliderCTD.sigma0;
    z = gliderCTD.sndspd;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    #xmin = NaNMath.mean(x) .- (nsig+0.5) * NaNMath.std(x);
    #xmax = NaNMath.mean(x) .+ (nsig+0.5) * NaNMath.std(x);
    #ymin = NaNMath.mean(y) .- (nsig+0.5) * NaNMath.std(y);
    #ymax = NaNMath.mean(y) .+ (nsig+0.5) * NaNMath.std(y);
    #zmin = NaNMath.mean(z) .- nsig * NaNMath.std(z);
    #zmax = NaNMath.mean(z) .+ nsig * NaNMath.std(z);

    xmin = pst.spice0min;
    xmax = pst.spice0max;
    ymin = pst.sigma0min;
    ymax = pst.sigma0max;
    zmin = pst.sndspdmin;
    zmax = pst.sndspdmax;

    fig = Figure(; size = tspres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Sigma0-Spice0",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Sigma0-Spice0",
        xlabel = "Spiciness",
        ylabel = "Potential Density Anomaly"
    )
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label="Sound Speed")
    fig
    #save(figoutdir * project * "_" * glidername * "_sigma0spice0.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_sigma0spice0.png", fig)


# plotting T/S diagram
    #x = saltAraw[1:pint:end]; 
    #y = ctempraw[1:pint:end];
    #z = sigma0raw[1:pint:end];
    x = gliderCTD.salt;
    y = gliderCTD.temp;
    z = gliderCTD.sigma0;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    #xmin = NaNMath.mean(x) .- (nsig+0.5) * NaNMath.std(x);
    #xmax = NaNMath.mean(x) .+ (nsig+0.5) * NaNMath.std(x);
    #ymin = NaNMath.mean(y) .- (nsig+0.5) * NaNMath.std(y);
    #ymax = NaNMath.mean(y) .+ (nsig+0.5) * NaNMath.std(y);
    #zmin = NaNMath.mean(z) .- nsig * NaNMath.std(z);
    #zmax = NaNMath.mean(z) .+ nsig * NaNMath.std(z);

    xmin = pst.saltmin;
    xmax = pst.saltmax;
    ymin = pst.tempmin;
    ymax = pst.tempmax;
    zmin = pst.sigma0min;
    zmax = pst.sigma0max;

    fig = Figure(; size = tspres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = project * " " * uppercasefirst(glidername) * " T/S",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " T/S",
        xlabel = "Practical Salinity",
        ylabel = "In-Situ Temperature"
    )
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label="Sigma0")
    fig
    #save(figoutdir * project * "_" * glidername * "_TSraw.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_TSraw.png", fig)

    #=
    # plotting Chl-a
    if isdefined(gliderCHLA, :t) == true
        x0 = datetime2unix.(DateTime("2022-01-01"));
        #x = chlatime;
        x = gliderCHLA.t;
        xdt, xtick, xticklabel = datetick(x);
        y = gliderCHLA.z;
        z = gliderCHLA.var;

        #findall(isnan.(electaCHLA.var)) 
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
    =#
end

function plotGliderCTD(gliderCTDarray, plotsetting, plotstruct)
    display("Plotting glider CTD data...")
    if (@isdefined plotsetting) == false
        ps = Glider.gliderPlotType.plotSetting[];
        for i = 1:length(gliderCTDarray)
            pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
            iday = 1; # day intervals for plotting
            ms = 6; # marker size
            tsms = 6; # time series marker size
            pres = (1600, 800); # plot resolution
            tspres = (1000, 1000); # time series plot resolution
            fs = 42; # font size
            ps = push!(ps, Glider.gliderPlotType.plotSetting(pint, iday, ms, tsms, pres, tspres, fs));
        end
    else
        ps = plotsetting;
    end

    if (@isdefined plotstruct) == false
        pst = Glider.gliderPlotType.plotStruct[];
        for i = 1:length(gliderCTDarray)
            figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/";
            project = "PASSENGERS"
            glidername = "Sylvia"
            latmin, latmax = 37, 40;
            lonmin, lonmax = -65.2, -59.8;            
            tempmin = 5.0;
            tempmax = 30.0;
            condmin = 3.0;
            condmax = 6.0;
            saltmin = 35.5;
            saltmax = 37.5;
            sigma0min = 20.0;
            sigma0max = 30.0;
            spice0min = -2.0;
            spice0max = 2.0;
            sndspdmin = 1450;
            sndspdmax = 1530;
            pst = push!(pst, Glider.gliderPlotType.plotStruct(figoutdir, project, glidername, tempmin, tempmax, condmin, condmax, saltmin, saltmax, sigma0min, sigma0max, spice0min, spice0max, sndspdmin, sndspdmax));
        end
    else
        pst = plotstruct;
    end

    for i = 1:length(gliderCTDarray)
        display(i)
        gliderCTDraw = gliderCTDarray[i];
        #lonrange = [NaNMath.minimum(gliderCTDraw.lon) NaNMath.maximum(gliderCTDraw.lon)];
        #latrange = [NaNMath.minimum(gliderCTDraw.lat) NaNMath.maximum(gliderCTDraw.lat)]; 

        #ctempstd = NaNMath.std(gliderCTDraw.ctemp);
        #condstd = NaNMath.std(gliderCTDraw.cond);
        #saltAstd = NaNMath.std(gliderCTDraw.saltA);
        #sigma0std = NaNMath.std(gliderCTDraw.sigma0);
        #sndspdstd = NaNMath.std(gliderCTDraw.sndspd);
        #spice0std = NaNMath.std(gliderCTDraw.spice0);

        #nsig = 2.5

        #figoutdir = "/Users/gong/GitHub/jlglider/slocum/figures/";
        #figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/";
        #ctemprange = (NaNMath.mean(gliderCTDraw.ctemp) .- nsig*ctempstd, NaNMath.mean(gliderCTDraw.ctemp) .+ nsig*ctempstd);
        #condrange = (NaNMath.mean(gliderCTDraw.cond) .- nsig*condstd, NaNMath.mean(gliderCTDraw.cond) .+ nsig*condstd);
        #saltArange = (NaNMath.mean(gliderCTDraw.saltA) .- nsig*saltAstd, NaNMath.mean(gliderCTDraw.saltA) .+ nsig*saltAstd);
        #sigma0range = (NaNMath.mean(gliderCTDraw.sigma0) .- nsig*sigma0std, NaNMath.mean(gliderCTDraw.sigma0) .+ nsig*sigma0std);
        #sndspdrange = (NaNMath.mean(gliderCTDraw.sndspd) .- nsig*sndspdstd, NaNMath.mean(gliderCTDraw.sndspd) .+ nsig*sndspdstd);
        #spice0range = (NaNMath.mean(gliderCTDraw.spice0) .- nsig*spice0std, NaNMath.mean(gliderCTDraw.spice0) .+ nsig*spice0std);

        #temprange = ctemprange;
        #saltrange = saltArange;
        #pst = plotStruct(figoutdir, gliderCTDraw.project, gliderCTDraw.glidername, temprange[1], temprange[2], condrange[1], condrange[2], saltrange[1], saltrange[2], sigma0range[1], sigma0range[2], spice0range[1], spice0range[2], sndspdrange[1], sndspdrange[2]);

        plot_glider_ctd(gliderCTDraw, ps[i], pst[i]);
    end
    display("Done.")
end

function plotGliderMap(gliderCTDarray, pst; pzrange=[-40, -30], varname="saltA", logzflag=0)
    #using NCDatasets, GLMakie, NaNMath, Statistics, Dates

    if (@isdefined logzflag) == false
        logzflag = 0;
    end

    bathypath = "/Users/gong/oceansensing Dropbox/C2PO/Data/bathy/ETOPO1/ETOPO_2022_v1_30s_N90W180_surface.nc";
    bathyds = Dataset(bathypath,"r");
    
    lon = bathyds["lon"][:];
    lat = bathyds["lat"][:];
    
    #latmin, latmax = 37, 40;
    #lonmin, lonmax = -65.2, -59.8;
    lonmin = pst[1].lonmin
    lonmax = pst[1].lonmax
    latmin = pst[1].latmin
    latmax = pst[1].latmax
    
    # approximate x-axis scaling to make it look "normal"
    dlat = latmax - latmin;
    dlon = lonmax - lonmin;
    lat0 = NaNMath.mean(gliderCTDarray[1].lat);
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
    
    plottitle = uppercase(gliderCTDarray[1].project) * " " * gliderCTDarray[1].glidertype * " " * varname * " (" * string(year(unix2datetime(NaNMath.minimum(gliderCTDarray[1].t)))) * "-" * string(year(unix2datetime(NaNMath.maximum(gliderCTDarray[1].t)))) * ")";
    plotname = uppercase(gliderCTDarray[1].project) * "_" * lowercase(gliderCTDarray[1].glidertype) * "_" * varname * ".png" 

    fig = Figure(size = pres, fontsize = 64)
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        xlabel = "Longitude",
        ylabel = "Latitude",
    )
    Makie.contourf!(x, y, zp, colormap = :bukavu, levels = range(zrange[1], zrange[2], length = 128))
    xlims!(lonmin, lonmax);
    ylims!(latmin, latmax);
    
    for ii = 1:length(gliderCTDarray)
        zind = findall(minimum(pzrange) .<= gliderCTDarray[ii].z .<= maximum(pzrange));
        local x = gliderCTDarray[ii].lon[zind];
        local y = gliderCTDarray[ii].lat[zind];
        local c = gliderCTDarray[ii].saltA[zind];
        cmin = minimum(c);
        cmax = maximum(c);
        display(ii)
        GLMakie.scatter!(x, y, color=c, colorrange=(cmin,cmax), colormap=:jet, markersize=20);
    end
    fig
    save(pst[1].figoutdir * plotname, fig)
    GLMakie.closeall()
end

end