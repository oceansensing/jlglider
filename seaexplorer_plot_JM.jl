#include("seaexplorer_load_rt.jl")
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt

using Plots, NaNMath, Dates
plotly()

figoutdir = "/Users/gong/Research/sea064/figures/"

# Jan Mayen
t1 = Dates.DateTime(2022,10,21,12,0,0);
t2 = Dates.DateTime(2022,10,31,12,0,0);

# Lofoten Basin Eddy
#t1 = Dates.DateTime(2022,11,03,12,0,0);
#t2 = Dates.DateTime(2022,11,30,12,0,0);

tind = findall(t1 .<= sea064pld1d.t .<= t2);

l8out6 = @layout([a; b; c; d; e; f])
l8out5 = @layout([a; b; c; d; e])
l8out4 = @layout([a; b; c; d])
l8out3 = @layout([a; b; c])

ms = 2;

htemp = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = temp[tind], seriestype=:scatter, c=:thermal, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.5, 6.0), colorbar = false)
#Plots.contour!(sea064pld1d.t, -sea064pld1d.z, sigma0);
hsalt = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = salt[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(34.6, 35.0), colorbar = true)
hsigma0 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = sigma0[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(27.4, 28.1), colorbar = false)
hUeast = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = ad2cp_Ueast[tind], seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.4, 0.4), colorbar = false)
hUnorth = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = ad2cp_Unorth[tind], seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.4, 0.4), colorbar = false)

hMReps1 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_eps1[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-12, -6))
hMReps2 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_eps2[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-12, -6))
hMRqc1 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = mr_qc1[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRqc2 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = mr_qc2[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRsh1std = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_sh1_std[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRsh2std = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = log10.(mr_sh2_std[tind]), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")

hchla = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = chla[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0.0, 0.6))
hbb700 = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = bb700[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0.0, 0.01))
hcdom = Plots.plot(sea064pld1d.t[tind], -sea064pld1d.z[tind], zcolor = cdom[tind], seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.5, 0.5))

norsephysplot = Plots.plot(htemp, hsalt, hsigma0, hMReps1, hUeast, hUnorth, layout = l8out6, size=(1500,1500), framestyle=:box, legend=:outertopright, title=["Temperature" "Salinity" "Sigma0" "TKE EPS1" "U (east)" "V (north)"]);
norseMRplot = Plots.plot(hMReps1, hMReps2, hMRsh1std, hMRsh2std, hMRqc1, hMRqc2, layout = l8out6, size=(1500,1500), framestyle=:box, legend=:outertopright, title=["TKE EPS1" "TKE EPS2" "Shear 1 STDDEV" "Shear 2 STDDEV" "QC1" "QC2"]);

#norseoptcplot = Plots.plot(htemp, hchla, hbb700, hcdom, layout = l8out4, size=(1500,1000), framestyle=:box, legend=:outertopright, title=["Temperature" "Chlorophyll-a" "BB 700" "CDOM"]);

Plots.savefig(norsephysplot, figoutdir * "norse_sea064_physics_JM.html");
Plots.savefig(norseMRplot, figoutdir * "norse_sea064_MR_JM.html");

#gr()
#norseEPS1plot = Plots.plot(hMReps1, size = (1000,800), framestyle=:box, markersize = 3, title="SEA064 RT EPS1")
#Plots.savefig(norseEPS1plot, figoutdir * "norse_sea064_MReps1.png");

Plots.savefig(norseoptcplot, figoutdir * "norse_sea064_optics_JM.html");

#norseplot = Plots.plot(htemp, hsalt, layout = l8out, size=(1300,1300), framestyle=:box, legend=:outertopright, title=["temperature" "salinity"]);
#Plots.savefig(norseplot, "norse_temp_salt.html")

#pFLUOR = Plots.plot(Dates.Date.(flt), flz, zcolor = flv, seriestype=:scatter, c=:algae, markersize = 3, markerstrokewidth = 0, legend = false, label="Chl-a Fluorescence");
#Plots.plot!(Dates.Date.([t1 t2]),[0 0], label="")
