#include("seaexplorer_load_rt.jl")
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt

using Plots, NaNMath, Dates, ColorSchemes, GibbsSeaWater
#plotly()
gr()

figoutdir = "/Users/gong/oceansensing Dropbox/C2PO/NORSE/sea064-20231112-norse/figures/"

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
elseif mission == 37 # Jan Mayen 2022
    region = "JM"
    #lims_temp = (-0.5, 6.0);
    lims_temp = (0.0, 8.0);
    lims_salt = (34.5, 35.2);
    lims_sigma0 = (27.4, 28.1);
    lims_spice0 = (-0.1, 0.55);
    lims_uv = (-0.6, 0.6);
    lims_sndspd = (1450, 1475);
    lims_n2 = (0, 0.00005);
    lims_epsilon = (-11, -5);
    t1 = Dates.DateTime(2022,10,21,12,0,0);
    t2 = Dates.DateTime(2022,10,31,12,0,0);
elseif mission == 48 # Jan Mayen 2023
    region = "JM"
    #lims_temp = (-0.5, 6.0);
    lims_temp = (-0.5, 4.5);
    lims_salt = (33.8, 35.2);
    lims_sigma0 = (26.9, 28.1);
    lims_spice0 = (-0.2, 0.5);
    lims_uv = (-0.6, 0.6);
    lims_sndspd = (1450, 1475);
    lims_n2 = (0, 0.00005);
    lims_epsilon = (-11.5, -5.5);
    t1 = Dates.DateTime(2023,11,12,0,0,0);
    t2 = Dates.DateTime(2023,11,27,0,0,0);
end

ms = 4;
ps = (1200,600);

#l8out8 = @layout([a; b; c; d; e; f; g; h])
#l8out7 = @layout([a; b; c; d; e; f; g])
#l8out6 = @layout([a; b; c; d; e; f])
#l8out5 = @layout([a; b; c; d; e])
#l8out4 = @layout([a; b; c; d])
#l8out3 = @layout([a; b; c])


saltAi = lims_salt[1]-0.5 : 0.005 : lims_salt[end]+0.5;
ctempi = lims_temp[1]-1 : 0.1 : lims_temp[end]+1;
ssaltA = repeat(reshape(saltAi, 1, :), length(ctempi), 1);
cctemp = repeat(ctempi, 1, length(saltAi));
ssigma0 = GibbsSeaWater.gsw_sigma0.(ssaltA, cctemp);
sspice0 = GibbsSeaWater.gsw_spiciness0.(ssaltA, cctemp);

revSpectral = reverse(ColorSchemes.Spectral);
loadcolorscheme(:harpercm, revSpectral.colors, revSpectral.notes);

tind = findall(t1 .<= sea064pld1d.t .<= t2);
tindmid = findall((t1 .<= sea064pld1d.tmid .<= t2) .&& (sea064pld1d.zmid .<= -2.0));

htemp = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.ctemp[tind], seriestype=:scatter, seriescolor=:harpercm, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_temp, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 Conservative Temperature");
#Plots.contour!(sea064pld1d.t, sea064pld1d.z, sigma0);
hsalt = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.saltA[tind], seriestype=:scatter, c=:haline, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_salt, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 Absolute Salinity")
hsigma0 = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.sigma0[tind], seriestype=:scatter, c=:dense, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_sigma0, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 Potential Density Anomaly")
hspice0 = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.spice0[tind], seriestype=:scatter, c=:hsv, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_spice0, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 Spiciness")
hsndspd = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.sndspd[tind], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_sndspd, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 Sound Speed")
hUeast = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.ad2cp_Ueast[tind], seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_uv, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 U: E-W Velocity")
hUnorth = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.ad2cp_Unorth[tind], seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_uv, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 V: N-S Velocity")

hN2 = Plots.plot(sea064pld1d.tmid[tindmid], sea064pld1d.zmid[tindmid], zcolor = sea064pld1d.n2[tindmid], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_n2, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 N^2")

eps1gind = findall(isnan.(log10.(sea064pld1d.mr_eps1[tind])) .== 0);
hMReps1 = Plots.plot(sea064pld1d.t[tind][eps1gind], sea064pld1d.z[tind][eps1gind], zcolor = log10.(sea064pld1d.mr_eps1[tind][eps1gind]), seriestype=:scatter, c=:harpercm, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_epsilon, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 TKE Dissipation: Shear 1")

hMReps2 = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = log10.(sea064pld1d.mr_eps2[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_epsilon, colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 TKE Dissipation: Shear 2")
hMRqc1 = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.mr_qc1[tind], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRqc2 = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.mr_qc2[tind], seriestype=:scatter, c=:jet1, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRsh1std = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = log10.(sea064pld1d.mr_sh1_std[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 Shear 1 std. dev.")
hMRsh2std = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = log10.(sea064pld1d.mr_sh2_std[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", colorbar = true, size=ps, dpi=300, framestyle=:box, title="SEA064 Shear 2 std. dev.")

#hchla = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.chla[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0, 0.6), colorbar = true, size=ps, framestyle=:box, title="Chlorophyll-a Fluorescence")
#hbb700 = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = log10.(sea064pld1d.bb700[tind] .+ 0.0001), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-5, -3), colorbar = true, size=ps, framestyle=:box, title="Backscattering (700 nm)")
#hcdom = Plots.plot(sea064pld1d.t[tind], sea064pld1d.z[tind], zcolor = sea064pld1d.cdom[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0, 1.0), colorbar = true, size=ps, framestyle=:box, title="CDOM Fluorescence")

hTS = Plots.plot(sea064pld1d.saltA[tind], sea064pld1d.ctemp[tind], zcolor = sea064pld1d.sigma0[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=lims_sigma0, xlims=lims_salt, ylims=lims_temp, colorbar = false, size=(1000,1000), framestyle=:box, title="SEA064 CT vs SA")
Plots.contour!(saltAi, ctempi, ssigma0, contour_labels=true, seriescolor = :black)
#Plots.contour!(saltAi, ctempi, sspice0, contour_labels=true, seriescolor = :black)

#norsephysplot = Plots.plot(htemp, hsalt, hsigma0, hspice0, hsndspd, hMReps1, hUeast, hUnorth, layout = l8out8, size=(1500,1700), framestyle=:box, legend=:outertopright, title=["Temperature" "Salinity" "Sigma0" "Spice0" "Sound Speed" "TKE EPS1" "U (east)" "V (north)"]);
#norseMRplot = Plots.plot(hN2, hMReps1, hMReps2, hMRsh1std, hMRsh2std, hMRqc1, hMRqc2, layout = l8out7, size=(1500,1500), framestyle=:box, legend=:outertopright, title=["N2" "TKE EPS1" "TKE EPS2" "Shear 1 STDDEV" "Shear 2 STDDEV" "QC1" "QC2"]);
#norseTSplot = Plots.plot(hTS, size=(1000,1000), framestyle=:box, title="CT vs SA");
#norseoptcplot = Plots.plot(htemp, hchla, hbb700, hcdom, layout = l8out4, size=(1500,1000), framestyle=:box, legend=:outertopright, title=["Temperature" "Chlorophyll-a" "BB 700" "CDOM"]);

if Plots.backend_name() != :plotly
    Plots.savefig(htemp, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_ctemp.png")
    Plots.savefig(hsalt, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_saltA.png")
    Plots.savefig(hsigma0, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_sigma0.png")
    Plots.savefig(hspice0, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_spice0.png")
    Plots.savefig(hsndspd, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_sndspd.png")
    Plots.savefig(hUeast, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_AD2CP_U.png")
    Plots.savefig(hUnorth, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_AD2CP_V.png")
    Plots.savefig(hTS, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_TS.png");
    Plots.savefig(hN2, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_N2.png")
    Plots.savefig(hMReps1, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_eps1.png")
    Plots.savefig(hMReps2, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_eps2.png")
    Plots.savefig(hMRsh1std, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_sh1stddev.png")
    Plots.savefig(hMRsh2std, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_MR1000G_sh2stddev.png")
    #Plots.savefig(hchla, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_chla.png")
    #Plots.savefig(h bb700, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_bb700.png")
    #Plots.savefig(hcdom, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_cdom.png")
    
    #Plots.savefig(norsephysplot, figoutdir * "norse_sea064_" * region * "_physics.png");
    #Plots.savefig(norseMRplot, figoutdir * "norse_sea064_" * region * "_MR.png");
else
    Plots.savefig(htemp, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_RBR_ctemp.html")
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
    #Plots.savefig(hchla, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_chla.html")
    #Plots.savefig(hbb700, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_bb700.html")
    #Plots.savefig(hcdom, figoutdir * "norse_sea064_M" * string(mission, pad=3) * "_" * region * "_FLBBCD_cdom.html")
    
    #Plots.savefig(norsephysplot, figoutdir * "norse_sea064_" * region * "_physics.html");
    #Plots.savefig(norseMRplot, figoutdir * "norse_sea064_" * region * "_MR.html");
end  