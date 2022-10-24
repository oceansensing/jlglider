include("seaexplorer_load_rt.jl")

using Plots, NaNMath, Dates
plotly()

t1 = Dates.DateTime(2022,10,21,12,0,0);
t2 = Dates.DateTime(2022,10,24,12,0,0);

function cleanEPS(epsin)
    epsout = deepcopy(epsin);
    badind = findall((epsin .>= 1e-4) .|| (epsin .<= 1e-13));
    epsout[badind] .= NaN;
    epsout = convert(Vector{Float64}, epsout)
    return epsout;
end

eps1 = cleanEPS(sea064pld1d.mr1000g_eps1);
eps2 = cleanEPS(sea064pld1d.mr1000g_eps2);

l8out = @layout([a; b])
heps1 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(eps1), seriestype=:scatter, c=:jet, markersize = 3, markerstrokewidth = 0, legend = false)
heps2 = Plots.plot(sea064pld1d.t, -sea064pld1d.z, zcolor = log10.(eps2), seriestype=:scatter, c=:jet, markersize = 3, markerstrokewidth = 0, legend = false)
norseplot = Plots.plot(heps1, heps2, layout = l8out, size=(1300,1300), framestyle=:box, legend=:outertopright, title=["SEA064 TKE EPS1" "SEA064 TKE EPS2"]);
Plots.savefig(norseplot, "norse_epsilon.html")

#pFLUOR = Plots.plot(Dates.Date.(flt), flz, zcolor = flv, seriestype=:scatter, c=:algae, markersize = 3, markerstrokewidth = 0, legend = false, label="Chl-a Fluorescence");
#Plots.plot!(Dates.Date.([t1 t2]),[0 0], label="")
