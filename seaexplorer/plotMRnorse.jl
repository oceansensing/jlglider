# this script plots the MR data collected during NORSE mission processed by Laur Ferris
# gong@vims.edu, 2024-11-13

using Distributed
using MAT
#using Plots
using GLMakie

include("/Users/gong/GitHub/ocean_julia/C2PO.jl")
import .C2PO: datetime2unix, unix2datetime, datenum2datetime, unix2yearday, yearday2datetime

if @isdefined(jm22) == false
    jm22path = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/mr1000g_processing/jm22.mat";
    jm22file = matopen(jm22path)
    jm22 = read(jm22file, "d")
    close(jm22file)
end

if @isdefined(jm23) == false
    jm23path = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20231112-norse-janmayen-complete/mr_processing/jm23.mat";
    jm23file = matopen(jm23path)
    jm23 = read(jm23file, "d")
    close(jm23file)
end

function plotEPS(MRdata, title, figoutpath, ylimits)
    t = datenum2datetime.(MRdata["time"][:]) 
    unixt = datetime2unix.(t);
    depth = -0.9935*MRdata["pres"][:];
    eps1 = MRdata["eps1"][:];

    pind = findall(unixt .> 1e5);
    t = t[pind];
    unixt = unixt[pind];
    depth = depth[pind];
    eps1 = eps1[pind];

    x = unix2yearday.(unixt);
    y = depth;
    c = eps1;

    cmin, cmax = -12, -6

    rmprocs(workers())
    fig = Figure(size = (1600, 800), fontsize = 28)
    ax = Axis(fig[1, 1],
        title = title,
        xlabel = "Year Day of " * string(year(t[1])),
        ylabel = "Depth (m)"
    )
    Makie.scatter!(ax, x, y, color = log10.(c), colormap = :jet, marker = :circle, markersize = 8, colorrange=(cmin, cmax))
    ylims!(ax, ylimits[1], ylimits[2])

    Colorbar(fig[1, 2], limits = (cmin, cmax), colormap = :jet, flipaxis = true, label = "Epsilon log(W/kg)")
    fig
    save(figoutpath, fig)
    GLMakie.closeall()
end

title22 = "TKE Dissipation (SEA064 NORSE 2022)";
figoutpath22 = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/ONR_annual_report_NORSE/2024/NORSE2022_MR.png"; 
plotEPS(jm22, title22, figoutpath22, [-1000, 0])

title22 = "TKE Dissipation (SEA064 NORSE 2022)";
figoutpath22 = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/ONR_annual_report_NORSE/2024/NORSE2022_MR_600.png"; 
plotEPS(jm22, title22, figoutpath22, [-600, 0])

title23 = "TKE Dissipation (SEA064 NORSE 2023)";
figoutpath23 = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/ONR_annual_report_NORSE/2024/NORSE2023_MR.png"; 
plotEPS(jm23, title23, figoutpath23, [-600, 0])

title23 = "TKE Dissipation (SEA064 NORSE 2023)";
figoutpath23 = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/ONR_annual_report_NORSE/2024/NORSE2023_MR_350.png"; 
plotEPS(jm23, title23, figoutpath23, [-350, 0])