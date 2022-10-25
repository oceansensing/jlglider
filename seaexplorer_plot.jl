#include("seaexplorer_load_rt.jl")
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt

using Plots, NaNMath, Dates
plotly()

t1 = Dates.DateTime(2022,10,21,12,0,0);
t2 = Dates.DateTime(2022,10,24,12,0,0);

l8out = @layout([a; b; c; d])

htemp = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = temp, seriestype=:scatter, c=:jet, markersize = 3, markerstrokewidth = 0, legend = false, label="")
hsalt = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = salt, seriestype=:scatter, c=:jet, markersize = 3, markerstrokewidth = 0, legend = false, label="")

heps1 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(eps1), seriestype=:scatter, c=:jet, markersize = 3, markerstrokewidth = 0, legend = false, label="")
heps2 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(eps2), seriestype=:scatter, c=:jet, markersize = 3, markerstrokewidth = 0, legend = false, label="")

norseplot = Plots.plot(htemp, hsalt, heps1, heps2, layout = l8out, size=(800,1000), framestyle=:box, legend=:outertopright, title=["Temperature" "Salinity" "SEA064 TKE EPS1" "SEA064 TKE EPS2"]);
Plots.savefig(norseplot, "norse_sea064.html")

#norseplot = Plots.plot(htemp, hsalt, layout = l8out, size=(1300,1300), framestyle=:box, legend=:outertopright, title=["temperature" "salinity"]);
#Plots.savefig(norseplot, "norse_temp_salt.html")

#pFLUOR = Plots.plot(Dates.Date.(flt), flz, zcolor = flv, seriestype=:scatter, c=:algae, markersize = 3, markerstrokewidth = 0, legend = false, label="Chl-a Fluorescence");
#Plots.plot!(Dates.Date.([t1 t2]),[0 0], label="")
