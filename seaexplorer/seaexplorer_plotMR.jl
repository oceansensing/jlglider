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

sh1file = ["jm22_eps1.csv"; "lbe22_eps1.csv"; "jm23_eps1.csv"];
sh2file = ["jm22_eps2.csv"; "lbe22_eps2.csv"; "jm23_eps2.csv"];
presfile = ["jm22_pressure.csv"; "lbe22_pressure.csv"; "jm23_pressure.csv"];
timefile = ["jm22_time.csv"; "lbe22_time.csv"; "jm23_time.csv"];

for i = 1:3
#i = 2
    display(i)
    if i == 1
        gdata = jm22;
    elseif i == 2
        gdata = lbe22;
    elseif i == 3
        gdata = jm23;
    end

    #gliderDataDir = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
    #gliderMissionDir = lowercase(glidername) * "-" * missiondate[i] * "-" * lowercase(projectname) * "-" * region[i] * "-complete";
    #gliderADCPdir = gliderDataDir * gliderMissionDir * "/ad2cp/" * "m" * missionid[i] * "_processed/";
    #gliderADCPpath = gliderADCPdir * "abs_dive_ocean_vel.csv";
    #figoutdir = gliderADCPdir * "figures/";

    gliderMRdir = "/Users/gong/oceansensing Dropbox/C2PO/NORSE/epsilon_csv/eps_csv/";
    gliderMReps1 = gliderMRdir * sh1file[i];
    gliderMReps2 = gliderMRdir * sh2file[i];
    gliderMRpres = gliderMRdir * presfile[i];
    gliderMRtime = gliderMRdir * timefile[i];

    #        df = CSV.read(navfilepath, header=1, delim=";", DataFrame, buffer_in_memory=true);

    eeps1 = Matrix(CSV.read(gliderMReps1, header=1, delim=",", DataFrame));
    eeps2 = Matrix(CSV.read(gliderMReps2, header=1, delim=",", DataFrame));
    ttstr = Matrix(CSV.read(gliderMRtime, header=1, delim=",", DataFrame));
    tt = Array{Float64, 2}(undef, size(ttstr, 1), size(ttstr, 2));
    pp = Matrix(CSV.read(gliderMRpres, header=1, delim=",", DataFrame));

    bdind = findall(eeps1 .== 0.0);
    eeps1[bdind] .= NaN;
    eeps2[bdind] .= NaN;
    pp[bdind] .= NaN;
    gdtind = findall(SubString.(ttstr,1,12) .!= "31-Dec--0001");
    bdtind = findall(SubString.(ttstr,1,12) .== "31-Dec--0001");
    tt[gdtind] = datetime2unix.(DateTime.(ttstr[gdtind], dateformat"d-u-y HH:MM:SS"));
    tt[bdtind] .= NaN;
    zz = -pp .* 0.993117;

    gind = findall(isnan.(tt) .== false);
    eps1 = eeps1[gind];
    eps2 = eeps2[gind];
    t = tt[gind];
    z = zz[gind];

#    yo = Int64.(gdf.yo_number);
#    t = collect(gdf.time_midpoint);
#    u = collect(gdf.u_ocean_vel);
#    v = collect(gdf.v_ocean_vel);
#    z = collect(gdf.depth_bins);

    #t0 = datetime2unix.(DateTime("2022-01-01"));
    #xdt, xtick, xticklabel = datetick(t);

    tdt = unix2datetime.(t); 
    yday = yearday.(tdt);

    pres = (1400, 600);
    fs = 28;
    ms = 10;
    zmin, zmax = -11, -5;

    figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/NORSE/OSM2024/";

    figeps1 = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figeps1[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " TKE Diss Shear 1 (post-processed)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yday, z, color=log10.(eps1), colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figeps1[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "log(W/kg)")
    figeps1
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_TKE_eps1.png", figeps1)
    GLMakie.closeall()

    figeps2 = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figeps2[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " TKE Diss Shear 2 (post-processed)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yday, z, color=log10.(eps2), colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figeps2[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "log(W/kg)")
    figeps2
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_TKE_eps2.png", figeps2)
    GLMakie.closeall()

    figeps1rt = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figeps1rt[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " TKE Diss Shear 1 (realtime)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yearday.(gdata.t), gdata.z, color=log10.(gdata.mr_eps1), colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figeps1rt[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "log(W/kg)")
    figeps1rt
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_TKE_eps1rt.png", figeps1rt)
    GLMakie.closeall()

    figeps2rt = Figure(resolution = pres, fontsize = fs)
    ax = Axis(figeps2rt[1, 1],
        title = projectname * " " * missiondate[i][1:4] * " " * glidername * " M" * missionid[i] * " TKE Diss Shear 2 (realtime)",
        xlabel = "Year Day",
        ylabel = "Depth (m)"
    )
    GLMakie.scatter!(yearday.(gdata.t), gdata.z, color=log10.(gdata.mr_eps2), colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
    Colorbar(figeps2rt[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = true, label = "log(W/kg)")
    figeps2rt
    save(figoutdir * projectname * "_" * glidername * "_M" * missionid[i] * "_" * region[i] * "_TKE_eps2rt.png", figeps2rt)
    GLMakie.closeall()
end