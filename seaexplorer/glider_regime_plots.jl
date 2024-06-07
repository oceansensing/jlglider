using Random
using Printf
using CairoMakie
using Dates
using GibbsSeaWater, Statistics, CSV, DataFrames

glider_data = CSV.read("/Users/rbbourdon/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20231112-norse-janmayen-complete/ad2cp/m48_processed/M48_science_data_processed.csv", DataFrame);
adcp_data = CSV.read("/Users/rbbourdon/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20231112-norse-janmayen-complete/ad2cp/m48_processed/absolute_ocean_vel.csv", DataFrame);

subset = glider_data[:, ["diveNum", "time", "depth", "SA", "CT"]]
subset1 = dropmissing(subset)
subset1[:, "time"] = subset1.time .* 10^-9
subset1[:, "TimeStamp"] = Dates.unix2datetime.(subset1.time);


wind_regime1 = DateTime("2023-11-12T00:00:00")
wind_regime2 = DateTime("2023-11-15T00:00:00")

bank1 = DateTime("2023-11-15T00:00:00")
bank2 = DateTime("2023-11-17T00:00:00")

wind_regime3 = DateTime("2023-11-17T00:00:00")
wind_regime4 = DateTime("2023-11-19T00:00:00")

bank3 = DateTime("2023-11-19T00:00:00")
bank4 = DateTime("2023-11-21T00:00:00")

storm_regime = DateTime("2023-11-21T00:00:00")
storm_regime2 = DateTime("2023-11-25T00:00:00")

jmch1 = DateTime("2023-11-25T00:00:00")
jmch2 = DateTime("2023-11-27T00:00:00")

regime1_df = filter(row -> wind_regime1 <= row.TimeStamp <= wind_regime2, subset1);
regime2_df = filter(row -> bank1 <= row.TimeStamp <= bank2, subset1);
regime3_df = filter(row -> wind_regime3 <= row.TimeStamp <= wind_regime4, subset1);
regime4_df = filter(row -> bank3 <= row.TimeStamp <= bank4, subset1);
storm_df = filter(row -> storm_regime <= row.TimeStamp <= storm_regime2, subset1);
jmch_df = filter(row -> jmch1 <= row.TimeStamp <= jmch2, subset1);

adcp_data[:, "TimeStamp"] = Dates.unix2datetime.(adcp_data.time_midpoint);

regime1_adcp = filter(row -> wind_regime1 <= row.TimeStamp <= wind_regime2, adcp_data);
regime2_adcp = filter(row -> bank1 <= row.TimeStamp <= bank2, adcp_data);
regime3_adcp = filter(row -> wind_regime3 <= row.TimeStamp <= wind_regime4, adcp_data);
regime4_adcp = filter(row -> bank3 <= row.TimeStamp <= bank4, adcp_data);
storm_adcp = filter(row -> storm_regime <= row.TimeStamp <= storm_regime2, adcp_data);
jmch_adcp = filter(row -> jmch1 <= row.TimeStamp <= jmch2, adcp_data);

ms = 3.5;
ps = (1000,500);
fs =14;

regime1_fig = Figure(resolution = ps, font = "Arial", fontsize = fs);
regime2_fig = Figure(resolution = ps, font = "Arial", fontsize = fs);
regime3_fig = Figure(resolution = ps, font = "Arial", fontsize = fs);
regime4_fig = Figure(resolution = ps, font = "Arial", fontsize = fs);
storm_fig = Figure(resolution = ps, font = "Arial", fontsize = fs);
jmch_fig = Figure(resolution = ps, font = "Arial", fontsize = fs);

ax_CT = Axis(regime1_fig[2, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 1 Conservative Temperature");
ax_SA = Axis(regime1_fig[2, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 1 Absolute Salinity");
ax_u = Axis(regime1_fig[3, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 1 Zonal Velocity");
ax_v = Axis(regime1_fig[3, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 1 Meridional Velocity");

ax_CT2 = Axis(regime2_fig[2, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 1 Conservative Temperature");
ax_SA2 = Axis(regime2_fig[2, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 1 Absolute Salinity");
ax_u2 = Axis(regime2_fig[3, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 1 Zonal Velocity");
ax_v2 = Axis(regime2_fig[3, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 1 Meridional Velocity");

ax_CT3 = Axis(regime3_fig[2, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 2 Conservative Temperature");
ax_SA3 = Axis(regime3_fig[2, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 2 Absolute Salinity");
ax_u3 = Axis(regime3_fig[3, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 2 Zonal Velocity");
ax_v3 = Axis(regime3_fig[3, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Wind Regime 2 Meridional Velocity");

ax_CT4 = Axis(regime4_fig[2, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 2 Conservative Temperature");
ax_SA4 = Axis(regime4_fig[2, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 2 Absolute Salinity");
ax_u4 = Axis(regime4_fig[3, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 2 Zonal Velocity");
ax_v4 = Axis(regime4_fig[3, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Bank 2 Meridional Velocity");

ax_CT5 = Axis(storm_fig[2, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Storm Conservative Temperature");
ax_SA5 = Axis(storm_fig[2, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Storm Absolute Salinity");
ax_u5 = Axis(storm_fig[3, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Storm Zonal Velocity");
ax_v5 = Axis(storm_fig[3, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Storm Meridional Velocity");

ax_CT6 = Axis(jmch_fig[2, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Jan Mayen Conservative Temperature");
ax_SA6 = Axis(jmch_fig[2, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Jan Mayen Absolute Salinity");
ax_u6 = Axis(jmch_fig[3, 1], xlabel = "Time", ylabel = "Depth (m)", title = "Jan Mayen Zonal Velocity");
ax_v6 = Axis(jmch_fig[3, 3], xlabel = "Time", ylabel = "Depth (m)", title = "Jan Mayen Meridional Velocity");

w1_CT = CairoMakie.scatter!(ax_CT, regime1_df.time, -regime1_df.depth, color=regime1_df.CT, colormap = :thermal);
Colorbar(regime1_fig[2,2], w1_CT; label = "Conservative Temperature (ºC)");
w1_SA = CairoMakie.scatter!(ax_SA, regime1_df.time, -regime1_df.depth, color=regime1_df.SA, markersize = ms, colormap = :viridis, colorrange = (34, 35.5));
Colorbar(regime1_fig[2,4], w1_SA; label = "Absolute Salinity (g/kg)");
w1_u = CairoMakie.scatter!(ax_u, regime1_adcp.time_midpoint, regime1_adcp.depth_bins, color=regime1_adcp.u_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime1_fig[3,2], w1_u; label = "Zonal Velocity (m/s)");
w1_v = CairoMakie.scatter!(ax_v, regime1_adcp.time_midpoint, regime1_adcp.depth_bins, color=regime1_adcp.v_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime1_fig[3,4], w1_v; label = "Meridional Velocity (m/s)");
regime1_timestamp = Date(2023, 11, 12):Day(1):Date(2023, 11, 15)
tick_positions = Float32.(LinRange(first(regime1_df.time), last(regime1_df.time), length(regime1_timestamp)))
tick_labels = string.(regime1_timestamp)
ax_CT.xticks = (tick_positions, tick_labels)
ax_SA.xticks = (tick_positions, tick_labels)
ax_u.xticks = (tick_positions, tick_labels)
ax_v.xticks = (tick_positions, tick_labels)
save("first_wind_regime_glider.png", regime1_fig)

b1_CT = CairoMakie.scatter!(ax_CT2, regime2_df.time, -regime2_df.depth, color=regime2_df.CT, markersize = ms, colormap = :thermal);
Colorbar(regime2_fig[2,2], b1_CT; label = "Conservative Temperature (ºC)");
b1_SA = CairoMakie.scatter!(ax_SA2, regime2_df.time, -regime2_df.depth, color=regime2_df.SA, markersize = ms, colormap = :viridis);
Colorbar(regime2_fig[2,4], b1_SA; label = "Absolute Salinity (g/kg)");
b1_u = CairoMakie.scatter!(ax_u2, regime2_adcp.time_midpoint, regime2_adcp.depth_bins, color=regime2_adcp.u_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime2_fig[3,2], b1_u; label = "Zonal Velocity (m/s)");
b1_v = CairoMakie.scatter!(ax_v2, regime2_adcp.time_midpoint, regime2_adcp.depth_bins, color=regime2_adcp.v_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime2_fig[3,4], b1_v; label = "Meridional Velocity (m/s)");
timestamp = Date(2023, 11, 15):Day(1):Date(2023, 11, 17)
tick_positions = Float32.(LinRange(first(regime2_df.time), last(regime2_df.time), length(timestamp)))
tick_labels = string.(timestamp)
ax_CT2.xticks = (tick_positions, tick_labels)
ax_SA2.xticks = (tick_positions, tick_labels)
ax_u2.xticks = (tick_positions, tick_labels)
ax_v2.xticks = (tick_positions, tick_labels)
save("first_bank_glider.png", regime2_fig)

w2_CT = CairoMakie.scatter!(ax_CT3, regime3_df.time, -regime3_df.depth, color=regime3_df.CT, markersize = ms, colormap = :thermal);
Colorbar(regime3_fig[2,2], w2_CT; label = "Conservative Temperature (ºC)");
w2_SA = CairoMakie.scatter!(ax_SA3, regime3_df.time, -regime3_df.depth, color=regime3_df.SA, markersize = ms, colormap = :viridis);
Colorbar(regime3_fig[2,4], w2_SA; label = "Absolute Salinity (g/kg)");
w2_u = CairoMakie.scatter!(ax_u3, regime3_adcp.time_midpoint, regime3_adcp.depth_bins, color=regime3_adcp.u_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime3_fig[3,2], w2_u; label = "Zonal Velocity (m/s)");
w2_v = CairoMakie.scatter!(ax_v3, regime3_adcp.time_midpoint, regime3_adcp.depth_bins, color=regime3_adcp.v_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime3_fig[3,4], w2_v; label = "Meridional Velocity (m/s)");
timestamp = Date(2023, 11, 17):Day(1):Date(2023, 11, 19)
tick_positions = Float32.(LinRange(first(regime3_df.time), last(regime3_df.time), length(timestamp)))
tick_labels = string.(timestamp)
ax_CT3.xticks = (tick_positions, tick_labels)
ax_SA3.xticks = (tick_positions, tick_labels)
ax_u3.xticks = (tick_positions, tick_labels)
ax_v3.xticks = (tick_positions, tick_labels)
save("second_wind_regime_glider.png", regime3_fig)

b2_CT = CairoMakie.scatter!(ax_CT4, regime4_df.time, -regime4_df.depth, color=regime4_df.CT, markersize = ms, colormap = :thermal);
Colorbar(regime4_fig[2,2], b2_CT; label = "Conservative Temperature (ºC)");
b2_SA = CairoMakie.scatter!(ax_SA4, regime4_df.time, -regime4_df.depth, color=regime4_df.SA, markersize = ms, colormap = :viridis);
Colorbar(regime4_fig[2,4], b2_SA; label = "Absolute Salinity (g/kg)");
b2_u = CairoMakie.scatter!(ax_u4, regime4_adcp.time_midpoint, regime4_adcp.depth_bins, color=regime4_adcp.u_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime4_fig[3,2], b2_u; label = "Zonal Velocity (m/s)");
b2_v = CairoMakie.scatter!(ax_v4, regime4_adcp.time_midpoint, regime4_adcp.depth_bins, color=regime4_adcp.v_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(regime4_fig[3,4], b2_v; label = "Meridional Velocity (m/s)");
timestamp = Date(2023, 11, 19):Day(1):Date(2023, 11, 21)
tick_positions = Float32.(LinRange(first(regime4_df.time), last(regime4_df.time), length(timestamp)))
tick_labels = string.(timestamp)
ax_CT4.xticks = (tick_positions, tick_labels)
ax_SA4.xticks = (tick_positions, tick_labels)
ax_u4.xticks = (tick_positions, tick_labels)
ax_v4.xticks = (tick_positions, tick_labels)
save("second_bank_glider.png", regime4_fig)

s_CT = CairoMakie.scatter!(ax_CT5, storm_df.time, -storm_df.depth, color=storm_df.CT, markersize = ms, colormap = :thermal);
Colorbar(storm_fig[2,2], s_CT; label = "Conservative Temperature (ºC)");
s_SA = CairoMakie.scatter!(ax_SA5, storm_df.time, -storm_df.depth, color=storm_df.SA, markersize = ms, colormap = :viridis, colorrange = (34, 35.5));
Colorbar(storm_fig[2,4], s_SA; label = "Absolute Salinity (g/kg)");
s_u = CairoMakie.scatter!(ax_u5, storm_adcp.time_midpoint, storm_adcp.depth_bins, color=storm_adcp.u_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(storm_fig[3,2], s_u; label = "Zonal Velocity (m/s)");
s_v = CairoMakie.scatter!(ax_v5, storm_adcp.time_midpoint, storm_adcp.depth_bins, color=storm_adcp.v_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(storm_fig[3,4], s_v; label = "Meridional Velocity (m/s)");
timestamp = Date(2023, 11, 21):Day(2):Date(2023, 11, 25)
tick_positions = Float32.(LinRange(first(storm_df.time), last(storm_df.time), length(timestamp)))
tick_labels = string.(timestamp)
ax_CT5.xticks = (tick_positions, tick_labels)
ax_SA5.xticks = (tick_positions, tick_labels)
ax_u5.xticks = (tick_positions, tick_labels)
ax_v5.xticks = (tick_positions, tick_labels)
save("storm_glider.png", storm_fig)

j_CT = CairoMakie.scatter!(ax_CT6, jmch_df.time, -jmch_df.depth, color=jmch_df.CT, markersize = ms, colormap = :thermal);
Colorbar(jmch_fig[2,2], j_CT; label = "Conservative Temperature (ºC)");
j_SA = CairoMakie.scatter!(ax_SA6, jmch_df.time, -jmch_df.depth, color=jmch_df.SA, markersize = ms, colormap = :viridis);
Colorbar(jmch_fig[2,4], j_SA; label = "Absolute Salinity (g/kg)");
j_u = CairoMakie.scatter!(ax_u6, jmch_adcp.time_midpoint, jmch_adcp.depth_bins, color=jmch_adcp.u_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(jmch_fig[3,2], j_u; label = "Zonal Velocity (m/s)");
j_v = CairoMakie.scatter!(ax_v6, jmch_adcp.time_midpoint, jmch_adcp.depth_bins, color=jmch_adcp.v_ocean_vel, markersize = ms, colormap = :balance, colorrange = (-0.5, 0.5));
Colorbar(jmch_fig[3,4], j_v; label = "Meridional Velocity (m/s)");
timestamp = Date(2023, 11, 25):Day(1):Date(2023, 11, 27)
tick_positions = Float32.(LinRange(first(jmch_df.time), last(jmch_df.time), length(timestamp)))
tick_labels = string.(timestamp)
ax_CT6.xticks = (tick_positions, tick_labels)
ax_SA6.xticks = (tick_positions, tick_labels)
ax_u6.xticks = (tick_positions, tick_labels)
ax_v6.xticks = (tick_positions, tick_labels)
save("jan_mayen_glider.png", jmch_fig)






