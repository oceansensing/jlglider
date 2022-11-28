#include("seaexplorer_load_rt.jl")
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt

using Plots, NaNMath, Dates, ColorSchemes
plotly()

figoutdir = "/Users/gong/Research/sea064/figures/"

if mission == 38 # LBE
    region = "LBE"
    #lims_temp = (-0.3, 8.7);
    lims_temp = (0.0, 8.0);
    lims_salt = (35.05, 35.41);
    lims_sigma0 = (27.2, 28.12);
    lims_spice0 = (-0.05, 1.05);
    lims_uv = (-0.6, 0.6);
    lims_sndspd = (1462, 1485);
    lims_n2 = (0, 0.00005);
    lims_epsilon = (-11, -5);
    t1 = Dates.DateTime(2022,11,02,12,0,0);
    t2 = Dates.DateTime(2022,11,30,12,0,0);
elseif mission == 37 # Jan Mayen
    region = "JM"
    #lims_temp = (-0.5, 6.0);
    lims_temp = (0.0, 8.0);
    lims_salt = (34.6, 35.2);
    lims_sigma0 = (27.4, 28.1);
    lims_spice0 = (-0.1, 0.55);
    lims_uv = (-0.6, 0.6);
    lims_sndspd = (1450, 1475);
    lims_n2 = (0, 0.00005);
    lims_epsilon = (-11, -5);
    t1 = Dates.DateTime(2022,10,21,12,0,0);
    t2 = Dates.DateTime(2022,10,31,12,0,0);
end

ms = 2.5;
ps = (1000,400);

tind = findall(t1 .<= sea064pld1d.t .<= t2);
tindmid = findall((t1 .<= tmid .<= t2) .&& (zmid .<= -2.0));

l8out8 = @layout([a; b; c; d; e; f; g; h])
l8out7 = @layout([a; b; c; d; e; f; g])
l8out6 = @layout([a; b; c; d; e; f])
l8out5 = @layout([a; b; c; d; e])
l8out4 = @layout([a; b; c; d])
l8out3 = @layout([a; b; c])


saltAi = lims_salt[1]-0.5 : 0.005 : lims_salt[end]+0.5;
ctempi = lims_temp[1]-1 : 0.1 : lims_temp[end]+1;
ssaltA = repeat(reshape(saltAi, 1, :), length(ctempi), 1);
cctemp = repeat(ctempi, 1, length(saltAi));
ssigma0 = gsw_sigma0.(ssaltA, cctemp);
sspice0 = gsw_spiciness0.(ssaltA, cctemp);

revSpectral = reverse(ColorSchemes.Spectral);
loadcolorscheme(:harpercm, revSpectral.colors, revSpectral.notes);

htemp = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = ctemp[tind], seriestype=:scatter, seriescolor=:harpercm, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_temp, colorbar = true, size=ps, dpi=300, framestyle=:box, title="Conservative Temperature");
#Plots.contour!(sea064pld1d.t, -sea064pld1d.z, sigma0);
hsalt = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = saltA[tind], seriestype=:scatter, c=:haline, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_salt, colorbar = true, size=ps, framestyle=:box, title="Absolute Salinity")
hsigma0 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = sigma0[tind], seriestype=:scatter, c=:dense, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_sigma0, colorbar = true, size=ps, framestyle=:box, title="Potential Density Anomaly")
hspice0 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = spice0[tind], seriestype=:scatter, c=:hsv, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_spice0, colorbar = true, size=ps, framestyle=:box, title="Spiciness")
hsndspd = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = sndspd[tind], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_sndspd, colorbar = true, size=ps, framestyle=:box, title="Sound Speed")
hUeast = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = ad2cp_Ueast[tind], seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_uv, colorbar = true, size=ps, framestyle=:box, title="U: E-W Velocity")
hUnorth = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = ad2cp_Unorth[tind], seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_uv, colorbar = true, size=ps, framestyle=:box, title="V: N-S Velocity")

hN2 = Plots.plot(tmid[tindmid], zmid[tindmid], zcolor = n2[tindmid], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_n2, colorbar = true, size=ps, framestyle=:box, title="N^2")
hMReps1 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_eps1[tind]), seriestype=:scatter, c=:harpercm, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_epsilon, colorbar = true, size=ps, framestyle=:box, title="TKE Dissipation: Shear 1")
hMReps2 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_eps2[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_epsilon, colorbar = true, size=ps, framestyle=:box, title="TKE Dissipation: Shear 2")
hMRqc1 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = mr_qc1[tind], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRqc2 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = mr_qc2[tind], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRsh1std = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_sh1_std[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", colorbar = true, size=ps, framestyle=:box, title="Shear 1 std. dev.")
hMRsh2std = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_sh2_std[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", colorbar = true, size=ps, framestyle=:box, title="Shear 2 std. dev.")

hchla = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = chla[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0, 0.6), colorbar = true, size=ps, framestyle=:box, title="Chlorophyll-a Fluorescence")
hbb700 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(bb700[tind] .+ 0.0001), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-5, -3), colorbar = true, size=ps, framestyle=:box, title="Backscattering (700 nm)")
hcdom = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = cdom[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0, 1.0), colorbar = true, size=ps, framestyle=:box, title="CDOM Fluorescence")

hTS = Plots.plot(saltA[tind], ctemp[tind], zcolor = sigma0[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_sigma0, xlims=lims_salt, ylims=lims_temp, colorbar = false, size=(1000,1000), framestyle=:box, title="CT vs SA")
Plots.contour!(saltAi, ctempi, ssigma0, contour_labels=true, seriescolor = :black)
#Plots.contour!(saltAi, ctempi, sspice0, contour_labels=true, seriescolor = :black)

#norsephysplot = Plots.plot(htemp, hsalt, hsigma0, hspice0, hsndspd, hMReps1, hUeast, hUnorth, layout = l8out8, size=(1500,1700), framestyle=:box, legend=:outertopright, title=["Temperature" "Salinity" "Sigma0" "Spice0" "Sound Speed" "TKE EPS1" "U (east)" "V (north)"]);
#norseMRplot = Plots.plot(hN2, hMReps1, hMReps2, hMRsh1std, hMRsh2std, hMRqc1, hMRqc2, layout = l8out7, size=(1500,1500), framestyle=:box, legend=:outertopright, title=["N2" "TKE EPS1" "TKE EPS2" "Shear 1 STDDEV" "Shear 2 STDDEV" "QC1" "QC2"]);
#norseTSplot = Plots.plot(hTS, size=(1000,1000), framestyle=:box, title="CT vs SA");
#norseoptcplot = Plots.plot(htemp, hchla, hbb700, hcdom, layout = l8out4, size=(1500,1000), framestyle=:box, legend=:outertopright, title=["Temperature" "Chlorophyll-a" "BB 700" "CDOM"]);

Plots.savefig(htemp, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_ctemp.html")
#Plots.savefig(htemp, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_ctemp.png")
Plots.savefig(hsalt, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_saltA.html")
Plots.savefig(hsigma0, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_sigma0.html")
Plots.savefig(hspice0, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_spice0.html")
Plots.savefig(hsndspd, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_sndspd.html")
Plots.savefig(hUeast, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_AD2CP_U.html")
Plots.savefig(hUnorth, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_AD2CP_V.html")
Plots.savefig(hTS, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_TS.html");
Plots.savefig(hN2, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_N2.html")
Plots.savefig(hMReps1, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_eps1.html")
Plots.savefig(hMReps2, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_eps2.html")
Plots.savefig(hMRsh1std, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_sh1stddev.html")
Plots.savefig(hMRsh2std, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_sh2stddev.html")
Plots.savefig(hchla, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_chla.html")
Plots.savefig(hbb700, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_bb700.html")
Plots.savefig(hcdom, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_cdom.html")

#Plots.savefig(norsephysplot, figoutdir * "norse_sea064_" * region * "_physics.html");
#Plots.savefig(norseMRplot, figoutdir * "norse_sea064_" * region * "_MR.html");
