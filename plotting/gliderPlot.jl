module gliderPlot

using Glider
include("/Users/gong/GitHub/jlglider/common/C2PO.jl")
using .C2PO: yearday2datetime, datetime2yearday
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
    ymin = pst.zmin;
    ymax = pst.zmax;
#    ymin = -350;
#    ymax = 0;

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

    #x = yday;
    x = tdt;
    #xdt, xtick, xticklabel = datetick(x);

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
        #xlabel = "Time",
        ylabel = "Depth (m)"
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:thermal, markersize=ms, colorrange=(zmin, zmax))
    #=
    ax.xticks = (xtick .- x0, xticklabel);
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#

    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :thermal, flipaxis = true, label = "Temperature (°C)")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_temp_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

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
        #xlabel = "Time",
        ylabel = "Depth (m)"
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Conductivity (S/m)")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_cond_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

    # plotting salinity
    y = gliderCTD.z;
    z = gliderCTD.salt;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = pst.saltmin;
    zmax = pst.saltmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Salinity",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Salinity",
        #xlabel = "Time",
        ylabel = "Depth (m)"
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:haline, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :haline, flipaxis = true, label = "Salinity (psu)")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_salt_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

    # plotting conservative temperature
    y = gliderCTD.z;
    z = gliderCTD.ctemp;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = pst.tempmin;
    zmax = pst.tempmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Conservative Temperature",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Conservative Temperature",
        #xlabel = "Time",
        ylabel = "Depth (m)"
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:thermal, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :thermal, flipaxis = true, label = "Conservative Temperature (°C)")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_ctemp_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

    # plotting absolute salinity
    #z = saltAraw[1:pint:end];
    y = gliderCTD.z;
    z = gliderCTD.saltA;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z); 
    zmin = pst.saltmin;
    zmax = pst.saltmax;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Absolute Salinity",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Absolute Salinity",
        #xlabel = "Time",
        ylabel = "Depth (m)"
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:haline, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :haline, flipaxis = true, label = "Absolute Salinity (g/kg)")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_saltA_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

    # plotting sigma0
    #z = sigma0raw[1:pint:end];
    y = gliderCTD.z;
    z = gliderCTD.sigma0;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    zmin = pst.sigma0min;
    zmax = pst.sigma0max;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = project * " " * uppercasefirst(glidername) * " Potential Density Anomaly",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Potential Density Anomaly",

        #xlabel = "Time",
        ylabel = "Depth (m)"
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:dense, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :dense, flipaxis = true, label = L"\sigma_{\theta} (kg/m^{3})")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_sigma0_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

    # plotting spice0
    #z = spice0raw[1:pint:end];
    y = gliderCTD.z;
    z = gliderCTD.spice0;
    #zmin = NaNMath.minimum(z);
    #zmax = NaNMath.maximum(z);
    zmin = pst.spice0min;
    zmax = pst.spice0max;
    fig = Figure(; size = pres, fontsize = fs)
    ax = Axis(fig[1, 1],
        #title = mission * " " * uppercasefirst(glidername) * " Spiciness0",
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Spice0",

        #xlabel = "Time",
        ylabel = "Depth (m)"
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #ax.xticks = (xtick .- x0, xticklabel);
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Spice0")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    #save(figoutdir * project * "_" * glidername * "_spice0.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_spice0_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

    # plotting turner angle
    SA = gliderCTD.saltA;
    CT = gliderCTD.ctemp; 
    p = gliderCTD.p*10; 
    Tu = Array{Float64, 1}(undef, length(SA)-1);
    Rsubrho = Array{Float64, 1}(undef, length(SA)-1); 
    p_mid = Array{Float64, 1}(undef, length(SA)-1);  
    GibbsSeaWater.gsw_turner_rsubrho(SA, CT, p, length(SA), Tu, Rsubrho, p_mid);
    c = Tu;
    cmin = -180;
    cmax = 180;
    y = GibbsSeaWater.gsw_z_from_p.(p_mid, NaNMath.mean(gliderCTD.lat), 0.0, 0.0); 

    jet8 = resample_cmap(:jet, 8);

    fig = Figure(; size = pres, fontsize = fs);
    ax = Axis(fig[1, 1],
        title = uppercase(project) * " " * yyyy0 * " " * uppercasefirst(glidername) * " Turner Angle",
        #xlabel = "Time",
        ylabel = "Depth (m)",
        #label = :bold
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x[1:end-1], y, color=c, colormap=jet8, markersize=10, colorrange=(cmin, cmax))
    cb = Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = jet8, flipaxis = true, label = "Turner Angle")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_Tu_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

    # plotting sound speed

    y = gliderCTD.z;
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
        #xlabel = "Time",
        ylabel = "Depth (m)",
        #label = :bold
    )
    ylims!(ax, ymin, ymax);
    Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    #=
    if length(xtick) > 10
        ax.xticks = (xtick[1:2:end] .- x0, xticklabel[1:2:end]);
    else
        ax.xticks = (xtick[1:1:end] .- x0, xticklabel[1:1:end]);
    end
    =#
    cb = Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "Sound Speed (m/s)")
    cb.labelrotation = 3*π/2 # Rotate by 90 degrees (π/2 radians)
    fig
    #save(figoutdir * project * "_" * glidername * "_soundspeed.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_soundspeed_" * string(Int(abs(ymin))) * ".png", fig)
    GLMakie.closeall()

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
    Makie.scatter!(x, y, color=z, colormap=:dense, markersize=tsms, colorrange=(zmin, zmax))
    #ax.xticks = (t[1]:86400:t[end], string.(Date.(td[1]:Day(1):td[end])))
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :dense, flipaxis = true, label="Sigma0")
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
    GLMakie.closeall()

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
    Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :dense, flipaxis = true, label="Sigma0")
    fig
    #save(figoutdir * project * "_" * glidername * "_TSraw.png", fig)
    save(figoutdir * project * "_" * glidername * "_" * yyyy0 * mm0 * dd0 * "_TSraw.png", fig)
    GLMakie.closeall()

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

function plotGliderTSarray(gliderCTDarray, ps, pst)
    display("Plotting plotGliderTSarray...")
    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/figures/";

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
    tmin = minimum(gliderCTDarray[1].t);
    tmax = maximum(gliderCTDarray[end].t);

    xmin = saltmin;
    xmax = saltmax;
    ymin = tempmin;
    ymax = tempmax;

    # time as color
    figTS = Figure(; size = tspres, fontsize = fs)
    ax = Axis(figTS[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " CT/SA - time",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature"
    )
    xmin = saltmin;
    xmax = saltmax;
    ymin = tempmin;
    ymax = tempmax;
    zmin = tmin;
    zmax = tmax;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.saltA;
        y = gliderCTD.ctemp;
        z = gliderCTD.t;
        Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    end
    Colorbar(figTS[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label="Unix Time")
    figTS
    save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_CA-SA-time.png", figTS);
    GLMakie.closeall()

    #density as color
    figTS = Figure(; size = tspres, fontsize = fs)
    ax = Axis(figTS[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " CT/SA - sigma0",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature"
    )
    xmin = saltmin;
    xmax = saltmax;
    ymin = tempmin;
    ymax = tempmax;
    zmin = sigma0min;
    zmax = sigma0max;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.saltA;
        y = gliderCTD.ctemp;
        z = gliderCTD.sigma0;
        Makie.scatter!(x, y, color=z, colormap=:dense, markersize=ms, colorrange=(zmin, zmax))
    end
    Colorbar(figTS[1, 2], limits = (zmin, zmax), colormap = :dense, flipaxis = true, label="sigma0 (kg/m^3)")
    figTS
    save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_CA-SA-sigma0.png", figTS);
    GLMakie.closeall()

    #sound speed as color
    figTS = Figure(; size = tspres, fontsize = fs)
    ax = Axis(figTS[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " CT/SA - c",
        xlabel = "Absolute Salinity",
        ylabel = "Conservative Temperature"
    )
    xmin = saltmin;
    xmax = saltmax;
    ymin = tempmin;
    ymax = tempmax;
    cmin = sndspdmin;
    cmax = sndspdmax;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.saltA;
        y = gliderCTD.ctemp;
        c = gliderCTD.sndspd;
        Makie.scatter!(x, y, color=c, colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
    end
    Colorbar(figTS[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="Sound Speed (m/s)")
    figTS
    save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_CA-SA-sndspd.png", figTS);
    GLMakie.closeall()

    display("Done plotting plotGliderTSarray.")

    #c/spice0, sigma0 as color 2D
    figSS = Figure(; size = tspres, fontsize = fs)
    ax = Axis(figSS[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " c/spice0 - sigma0",
        xlabel = "Spice0",
        ylabel = "Sounds Speed (m/s)"
    )
    xmin = spice0min;
    xmax = spice0max;
    ymin = 1470;
    ymax = 1550;
    cmin = sigma0min;
    cmax = sigma0max;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.spice0;
        y = gliderCTD.sndspd;
        c = gliderCTD.sigma0;
        Makie.scatter!(x, y, color=c, colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
    end
    Colorbar(figSS[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="Sound Speed (m/s)")
    figSS
    save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_SS-Spice0-sigma0.png", figSS);
    GLMakie.closeall()

    #c/spice0/sigma0 as color 3D
    figSS3d = Figure(; size = tspres, fontsize = fs)
    ax = Axis3(figSS3d[1, 1],
        title = uppercase(project) * " " * yyyy0 * "-" * yyyyN * " VIMS Gliders" * " c/spice0/sigma0 - temp",
        xlabel = "Spice0",
        ylabel = "Sounds Speed (m/s)",
        zlabel = "Sigma0 (kg/m^3)"
    )
    xmin = spice0min;
    xmax = spice0max;
    ymin = sndspdmin;
    ymax = sndspdmax;
    zmin = sigma0min;
    zmax = sigma0max;
    cmin = tempmin;
    cmax = tempmax;
    xlims!(ax, xmin, xmax);
    ylims!(ax, ymin, ymax);
    zlims!(ax, zmin, zmax);

    for i = 1:length(gliderCTDarray)
        gliderCTD = gliderCTDarray[i];
        x = gliderCTD.spice0;
        y = gliderCTD.sndspd;
        z = gliderCTD.sigma0;
        c = gliderCTD.ctemp;
        Makie.scatter!(x, y, z, color=c, colormap=:jet, markersize=ms, colorrange=(cmin, cmax))
    end
    Colorbar(figSS3d[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label="Conservative Temperature")
    figSS3d
    #save(figoutdir * project * "_glider_" * yyyy0 * "-" * yyyyN * "_SS-Spice0-sigma0-ctemp.png", figSS3d);
    #GLMakie.closeall()
    display("Done plotting plotGliderTSarray.")

end

function plotGliderMap(gliderCTDarray, pst; pzrange=[-30, -20], varname="saltA", logzflag=0)
    display("Plotting plotGliderMap...")
    #using NCDatasets, GLMakie, NaNMath, Statistics, Dates

    if (@isdefined logzflag) == false
        logzflag = 0;
    end

    i = length(gliderCTDarray);

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
        log10zb[pzbind] .= log10.(zb[pzbind]);
        log10zb[nzbind] .= -log10.(-zb[nzbind]);
        log10zb[zzbind] .= 0;
        zrange = (-4, 4);
        zbp = log10zb;
    else
        zbp = zb;
        zrange = (-6000, 6000);
    end
    

    plottitle = uppercase(gliderCTDarray[i].project) * " " * gliderCTDarray[i].glidertype * " " * varname * " (" * string(year(unix2datetime(NaNMath.minimum(gliderCTDarray[1].t)))) * "-" * string(year(unix2datetime(NaNMath.maximum(gliderCTDarray[end].t)))) * ")";
    plotname = uppercase(gliderCTDarray[i].project) * "_" * lowercase(gliderCTDarray[i].glidertype) * "_" * varname * ".png" 

    fig = Figure(size = pres, fontsize = 64)
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        xlabel = "Longitude",
        ylabel = "Latitude",
    )
    Makie.contourf!(xb, yb, zbp, colormap = ColorSchemes.bukavu, levels = range(zrange[1], zrange[2], length = 128))
    xlims!(lonmin, lonmax);
    ylims!(latmin, latmax);
    
    for ii = 1:length(gliderCTDarray)
        zind = findall(minimum(pzrange) .<= gliderCTDarray[ii].z .<= maximum(pzrange));
        local x = gliderCTDarray[ii].lon[zind];
        local y = gliderCTDarray[ii].lat[zind];
        local z = gliderCTDarray[ii].z[zind];

        if varname == "saltA"
            local c = gliderCTDarray[ii].saltA[zind];
            cmin = pst[ii].saltmin;
            cmax = pst[ii].saltmax;
            cmap = ColorSchemes.haline;
        elseif varname == "ctemp"
            local c = gliderCTDarray[ii].ctemp[zind];
            cmin = pst[ii].tempmin;
            cmax = pst[ii].tempmax;
            cmap = ColorSchemes.thermal;
        elseif varname == "spice0"
            local c = gliderCTDarray[ii].spice0[zind];
            cmin = pst[ii].spice0min;
            cmax = pst[ii].spice0max;
            cmap = ColorSchemes.jet;
        end
        display(ii)
        GLMakie.scatter!(x, y, color=c, colorrange=(cmin,cmax), colormap=cmap, markersize=20);
        #Colorbar(fig[1, 2], limits = (cmin,cmax), colormap=cmap, flipaxis=true, label=varname)
    end
    fig
    save(pst[i].figoutdir * plotname, fig)
    GLMakie.closeall()
    display("Done plotGliderMap.")
end

end