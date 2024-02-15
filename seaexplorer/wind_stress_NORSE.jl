using Glob, DataFrames, CSV, Dates, Missings, NaNMath, Interpolations
using NCDatasets
using GibbsSeaWater
using Plots, DSP, MAT, Statistics
plotly()

include("/Users/rbbourdon/Github/jlglider/seaexplorer/C2PO.jl")
include("/Users/rbbourdon/Github/jlglider/seaexplorer/as_consts.jl")
using .C2PO: stresslp

norse_era = NCDataset("/Users/rbbourdon/oceansensing Dropbox/C2PO/Data/ERA5/NORSE/2022-2023_10_11_12.nc", "r");
jm23_gli = CSV.read("/Users/rbbourdon/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20231112-norse-janmayen-complete/ad2cp/m48_processed/M48_glider_data_processed.csv", DataFrame);
jm22_gli = CSV.read("/Users/rbbourdon/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/ad2cp/m37_processed/M37_glider_data_processed.csv", DataFrame);
lbe_gli = CSV.read("/Users/rbbourdon/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/ad2cp/m38_processed/M38_glider_data_processed.csv", DataFrame);

u10 = norse_era["u10"];
v10 = norse_era["v10"];
lat = norse_era["latitude"][:, 1];
lon = norse_era["longitude"][:, 1];
time = norse_era["time"][:, 1];
jm23_gli.date_time = Dates.unix2datetime.(jm23_gli.date_float);
jm22_gli.date_time = Dates.unix2datetime.(jm22_gli.date_float);
lbe_gli.date_time = Dates.unix2datetime.(lbe_gli.date_float);


function find_nearest_value(arr, value)
    idx = argmin(abs.(arr .- value))
    return idx
end
using Dates

function find_nearest_time(arr, value)
    idx = argmin(abs.(arr .- value))
    return idx
end

function tau_timeseries(data::DataFrame, u10, v10, lat::Array, lon::Array, time::Array)

    gdf = groupby(data, :diveNum);
    dives = unique(data.diveNum);
    time_array = [];
    tau_array = [];
    total_wind_speed = [];
    for dive in dives
        yo = gdf[dive];
        start_row = yo[1, :];
        glider_lat = start_row.Lat_dd;
        glider_lon = start_row.Lon_dd;
        glider_time = start_row.date_time;
        era_lat_idx = find_nearest_value(lat, glider_lat);
        era_lon_idx = find_nearest_value(lon, glider_lon);
        era_time_idx = find_nearest_time(time, glider_time);
        yo_uwind = u10[era_lon_idx, era_lat_idx, 1, era_time_idx];
        yo_vwind = v10[era_lon_idx, era_lat_idx, 1, era_time_idx];
        tot_wind = sqrt(yo_uwind^2 + yo_vwind^2);
        z = 10;
        rhoa = 1.22;
        tau = stresslp(tot_wind, z, rhoa);
        push!(time_array, glider_time);
        push!(tau_array, tau);
        push!(total_wind_speed, tot_wind);
    end 
    return time_array, tau_array, total_wind_speed
end

jm23_time_array, jm23_tau, jm23_wind = tau_timeseries(jm23_gli, u10, v10, lat, lon, time);
jm22_time_array, jm22_tau, jm22_wind = tau_timeseries(jm22_gli, u10, v10, lat, lon, time);
lbe_time_array, lbe_tau, lbe_wind = tau_timeseries(lbe_gli, u10, v10, lat, lon, time);

jm23_time_array = DateTime.(jm23_time_array);
jm22_time_array = DateTime.(jm22_time_array);
lbe_time_array = DateTime.(lbe_time_array);

jm23_tau = [x[1] for x in jm23_tau];
jm22_tau = [x[1] for x in jm22_tau];
lbe_tau = [x[1] for x in lbe_tau];

jm23_wind = [x[1] for x in jm23_wind];
jm22_wind = [x[1] for x in jm22_wind];
lbe_wind = [x[1] for x in lbe_wind];



jm23_tau_timeseries = Plots.plot(jm23_time_array, jm23_tau, xlabel="Time", ylabel="Total Wind Stress (N/m^2)", title="Time Series of Wind Stress for NORSE M48", lw=2, grid=:on, size=(1100, 250), legend=false, framestyle=:box);
jm22_tau_timeseries = Plots.plot(jm22_time_array, jm22_tau, xlabel="Time", ylabel="Total Wind Stress (N/m^2)", title="Time Series of Wind Stress for NORSE M37", lw=2, grid=:on, size=(900, 250), legend=false, framestyle=:box);
lbe_tau_timeseries = Plots.plot(lbe_time_array, lbe_tau, xlabel="Time", ylabel="Total Wind Stress (N/m^2)", title="Time Series of Wind Stress for NORSE M38", lw=2, grid=:on, size=(900, 250), legend=false, framestyle=:box);

figoutdir = "/Users/rbbourdon/Desktop/NORSE/science_figs/"
Plots.savefig(jm23_tau_timeseries, figoutdir * "jm23_tau_timeseries.html");
Plots.savefig(jm22_tau_timeseries, figoutdir * "jm22_tau_timeseries.html");
Plots.savefig(lbe_tau_timeseries, figoutdir * "lbe_tau_timeseries.html")


