#include("seaexplorer_load_rt.jl")
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt

using Plots, NaNMath, Dates
plotly()

figoutdir = "/Users/gong/Research/sea064/figures/"

t1 = Dates.DateTime(2022,10,21,12,0,0);
t2 = Dates.DateTime(2022,10,24,12,0,0);

l8out6 = @layout([a; b; c; d; e; f])
l8out5 = @layout([a; b; c; d; e])
l8out4 = @layout([a; b; c; d])
l8out3 = @layout([a; b; c])

ms = 2;

htemp = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = temp, seriestype=:scatter, c=:thermal, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.5, 6), colorbar = false)
#Plots.contour!(sea064pld1d.t, -sea064pld1d.z, sigma0);
hsalt = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = salt, seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(34.6, 35), colorbar = true)
hsigma0 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = sigma0, seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(27.4, 28.1), colorbar = false)
hUeast = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = ad2cp_Ueast, seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.4, 0.4), colorbar = false)
hUnorth = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = ad2cp_Unorth, seriestype=:scatter, c=:vik, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.4, 0.4), colorbar = false)

hMReps1 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(mr_eps1), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMReps2 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(mr_eps2), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRqc1 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = mr_qc1, seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRqc2 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = mr_qc2, seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRsh1std = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(mr_sh1_std), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")
hMRsh2std = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(mr_sh2_std), seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="")

hchla = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = chla, seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0.0, 0.6))
hbb700 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = bb700, seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(0.0, 0.01))
hcdom = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = cdom, seriestype=:scatter, c=:jet, markersize = ms, markerstrokewidth = 0, legend = false, label="", clims=(-0.5, 0.5))

norsephysplot = Plots.plot(htemp, hsalt, hsigma0, hMReps1, hUeast, hUnorth, layout = l8out6, size=(1500,1500), framestyle=:box, legend=:outertopright, title=["Temperature" "Salinity" "Sigma0" "TKE EPS1" "U (east)" "V (north)"]);
norseMRplot = Plots.plot(hMReps1, hMReps2, hMRqc1, hMRqc2, hMRsh1std, hMRsh2std, layout = l8out6, size=(1500,1500), framestyle=:box, legend=:outertopright, title=["TKE EPS1" "TKE EPS2" "QC1" "QC2" "Shear 1 STDDEV" "Shear 2 STDDEV"]);

#norseoptcplot = Plots.plot(htemp, hchla, hbb700, hcdom, layout = l8out4, size=(1500,1000), framestyle=:box, legend=:outertopright, title=["Temperature" "Chlorophyll-a" "BB 700" "CDOM"]);

Plots.savefig(norsephysplot, figoutdir * "norse_sea064_physics.html");
Plots.savefig(norseMRplot, figoutdir * "norse_sea064_MR.html");

gr()
norseEPS1plot = Plots.plot(hMReps1, size = (1000,800), framestyle=:box, markersize = 3, title="SEA064 RT EPS1")
Plots.savefig(norseEPS1plot, figoutdir * "norse_sea064_MReps1.png");

#Plots.savefig(norseoptcplot, figoutdir * "norse_sea064_optics.html");

#norseplot = Plots.plot(htemp, hsalt, layout = l8out, size=(1300,1300), framestyle=:box, legend=:outertopright, title=["temperature" "salinity"]);
#Plots.savefig(norseplot, "norse_temp_salt.html")

#pFLUOR = Plots.plot(Dates.Date.(flt), flz, zcolor = flv, seriestype=:scatter, c=:algae, markersize = 3, markerstrokewidth = 0, legend = false, label="Chl-a Fluorescence");
#Plots.plot!(Dates.Date.([t1 t2]),[0 0], label="")
