using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations
using GLMakie, ColorSchemes

#GLMakie.activate!()

import slocumFunc: datetick
import slocumType: plotSetting, ctdStruct, sciStruct
import slocumFunc: pyrow2jlcol, intersectalajulia2, glider_ctd_load, glider_var_load, glider_presfunc

dbdreader = pyimport("dbdreader");
gsw = GibbsSeaWater;

datamode = "realtime" # delayed or realtime
mission = "PASSENGERS 2023";

#dataroot = "/mnt/c/Users/C2PO/oceansensing Dropbox/C2PO/";
dataroot = "/Users/gong/oceansensing Dropbox/C2PO/";

glidername_electa = "electa";
rootdir_electa = dataroot * "PASSENGERS/2023_glider_data/electa-20230523-passengers/";
fromgliderdir_electa = rootdir_electa * "from-glider/"; 
datadir_electa = fromgliderdir_electa * datamode * "/" * "electa-from-glider-20230612T200456/";
cacdir_electa = fromgliderdir_electa * "cache/";
figoutdir_electa = rootdir_electa * "figures/";

glidername = glidername_electa;
rootdir = rootdir_electa;
fromgliderdir = fromgliderdir_electa; 
#datadirpath = Glob.glob("electa-from-glider*.zip", fromgliderdir)[1];
datadir = datadir_electa;
cacdir = cacdir_electa;
figoutdir = figoutdir_electa;

pint = 1; # this is the data decimation for plotting. Makie is so fast that it's not necessary, but Plots.jl would need it. Not using Plots.jl because of a bug there with colormap
iday = 1; # day intervals for plotting
ms = 6;
tsms = 5;
pres = (1200, 800)
tspres = (1000, 1000)
ps = plotSetting(pint, iday, ms, tsms, pres, tspres);

# specify valid data time period
t0 = DateTime("2023-05-23");
tN = DateTime("2023-06-20");
trange = datetime2unix.([t0; tN]);

# setup glider data loading using dbdreader
if datamode == "realtime"
    dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files_only = true, cacheDir = cacdir, return_nan=true);
    dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", cacheDir = cacdir);
else
    dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files_only = true, cacheDir = cacdir, return_nan=true);
end
engvars = dataGlider.parameterNames["eng"];
scivars = dataGlider.parameterNames["sci"];

t, sci_m_present_time, lon, lat, m_pressure, sci_water_pressure, sci_water_temp, sci_water_cond  = dataGlider.get_sync("sci_m_present_time", "m_lon", "m_lat", "m_pressure", "sci_water_pressure", "sci_water_temp", "sci_water_cond");
t2, sci_m_present_time2, sci_flbbcd_chlor_units, sci_flbbcd_cdom_units, sci_flbbcd_bb_units, sci_bsipar_par =  dataGlider.get_sync("sci_m_present_time", "sci_flbbcd_chlor_units", "sci_flbbcd_cdom_units", "sci_flbbcd_bb_units", "sci_bsipar_par");
#tisorted = sortperm(sci_m_present_time);

mlon = NaNMath.mean(lon);
mlat = NaNMath.mean(lat);

t1, p1, temp1, cond1 = glider_ctd_load(sci_m_present_time, sci_water_pressure, sci_water_temp, sci_water_cond, trange);
t2, p2, temp2, temp2ind = glider_var_load(sci_m_present_time, sci_water_pressure, sci_water_temp, trange, [0.1 40.0]);
z1 = gsw.gsw_z_from_p.(p1*10, mlat, 0.0, 0.0);  
z2 = gsw.gsw_z_from_p.(p2*10, mlat, 0.0, 0.0)

#t, sci_m_present_time, m_pressure, sci_water_presure, temp_sci = dataGlider.get_sync("sci_m_present_time", "m_pressure", "sci_water_pressure", "sci_water_temp");
#tisorted = sortperm(t_sci);
#x = t_sci[tisorted][end-10000:end-5000];
#y = p_sci[tisorted][end-10000:end-5000];
#z = temp_sci[tisorted][end-10000:end-5000];

x = t1;
y = z1;
z = temp1;

xdt, xtick, xticklabel = datetick(x);

#GLMakie.plot(x, y, color=z)
x0 = datetime2unix.(DateTime("2022-01-01"))
#GLMakie.plot(x .- x0, y, color=z)

# plotting conservative temperature
zmin = 18.0;
zmax = 24.0; 
fig = Figure(resolution = pres)
ax = Axis(fig[1, 1],
    title = mission * " " * glidername * " Temperature",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
ax.xticks = (xtick .- x0, xticklabel);
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
fig
save(figoutdir * mission * "_" * glidername * "_temp1.png", fig)

#=
# plotting chlorophyll-a
z = sci_flbbcd_chlor_units;
zmin = 0;
zmax = 1.0; 
fig = Figure(resolution = pres)
ax = Axis(fig[1, 1],
    title = mission * " " * glidername * " FLBBCD Chl-a",
    xlabel = "Time",
    ylabel = "Depth"
)
Makie.scatter!(x .- x0, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
ax.xticks = (xtick, xticklabel);
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
fig
save(figoutdir * mission * "_" * glidername * "_chla.png", fig)

#Plots.scatter(x, y, zcolor=z, clim=(18,20), markerstrokewidth = 0)

#=
# load engineering data from raw glider DBD files
m_present_time = dataGlider.get("m_present_time");
m_lat = dataGlider.get("m_lat"); 
m_lon = dataGlider.get("m_lon"); 
m_pressure= dataGlider.get("m_pressure");
m_vacuum = dataGlider.get("m_vacuum");
m_battery = dataGlider.get("m_battery");
#m_leak_detect = dataGlider.get("m_leakdetect"); 
m_veh_temp = dataGlider.get("m_veh_temp");  
m_roll = dataGlider.get("m_roll");
c_heading = dataGlider.get("c_heading"); 
m_heading = dataGlider.get("m_heading"); 
m_pitch = dataGlider.get("m_pitch");
m_battpos = dataGlider.get("m_battpos");
m_de_oil_vol = dataGlider.get("m_de_oil_vol");
m_ballast_pumped = dataGlider.get("m_ballast_pumped");
m_fin = dataGlider.get("m_fin");
m_depth_rate = dataGlider.get("m_depth_rate"); 
m_altimeter_status = dataGlider.get("m_altimeter_status"); 
m_raw_altitude = dataGlider.get("m_raw_altitude");
m_altimeter_voltage = dataGlider.get("m_altimeter_voltage"); 
m_num_tot_inflections = dataGlider.get("m_tot_num_inflections"); 

#m_ = pyrow2jlcol(dataGlider.get("m_")); 

# load CTD data from raw glider EBD files
#tctd, cond, temp, pres, m_de_oil_vol = dataGlider.get_CTD_sync("m_de_oil_vol");
sci_m_present_time = dataGlider.get("sci_m_present_time");
sci_water_temp = dataGlider.get("sci_water_temp");
sci_water_cond = dataGlider.get("sci_water_cond");
sci_water_pressure = dataGlider.get("sci_water_pressure");
sci_ctd41cp_timestamp = dataGlider.get("sci_ctd41cp_timestamp"); 
ctdtime = sci_ctd41cp_timestamp[2];


sci_flbbcd_chlor_units = dataGlider.get("sci_flbbcd_chlor_units");
sci_flbbcd_cdom_units = dataGlider.get("sci_flbbcd_cdom_units");
sci_flbbcd_bb_units = dataGlider.get("sci_flbbcd_bb_units");
sci_bsipar_par = dataGlider.get("sci_bsipar_par");

# calculate derived values from CTD data
llat = Statistics.mean(m_lat[2]);
llon = Statistics.mean(m_lon[2]);
#llon = -73.4;
#llat = 38.0;

presfunc, presraw = glider_presfunc(sci_water_pressure, trange);
tempraw, temptime, temppres, tempz = glider_var_load(sci_water_temp, trange, [0.1 40.0], sci_water_pressure, llat)
condraw, condtime, condpres, condz = glider_var_load(sci_water_cond, trange, [0.01 100.0], sci_water_pressure, llat)

if isempty(sci_flbbcd_chlor_units) != true
    chlaraw, chlatime, chlapres, chlaz = glider_var_load(sci_flbbcd_chlor_units, trange, [-0.1 3.0], sci_water_pressure, llat)
end

if isempty(sci_flbbcd_cdom_units) != true
    cdomraw, cdomtime, cdompres, cdomz = glider_var_load(sci_flbbcd_cdom_units, trange, [-5.0 5.0], sci_water_pressure, llat)
end

if isempty(sci_flbbcd_bb_units) != true
    bb700raw, bb700time, bb700pres, bb700z = glider_var_load(sci_flbbcd_bb_units, trange, [0.0 0.008], sci_water_pressure, llat)
end

if isempty(sci_bsipar_par) != true
    bparraw, bpartime, bparpres, bparz = glider_var_load(sci_bsipar_par, trange, [0.0 6000.0], sci_water_pressure, llat)
end

# find common glider values
tctd = unique(intersect(presraw[:,1], tempraw[:,1], condraw[:,1]));
tctdT = intersectalajulia2(tctd, tempraw[:,1])[3];
tctdC = intersectalajulia2(tctd, condraw[:,1])[3];
tctdP = intersectalajulia2(tctd, presraw[:,1])[3];

#=
tpuck = chlaraw[:,1];

# build a common timeline (not necessary if data is of good quality and CTD data all has the same length)
tall = sort(unique(intersect(tempraw[:,1], condraw[:,1], presraw[:,1])));
tind = findall(trange[1] .<= tall .<= trange[end]);

t1 = floor(NaNMath.minimum(tall[tind]));
t2 = ceil(NaNMath.maximum(tall[tind]));
tt = collect(t1:1.0:t2);

# extract common time values for all CTD parameters using interpolated to a regularly spaced time grid
dtctd_fit = unix2datetime.(tall[tind]);
temp_fit = tempfunc(tall[tind]);
cond_fit = condfunc(tall[tind]);
pres_fit = presfunc(tall[tind]);
chla_fit = chlafunc(tall[tind]);
tctd_fit = deepcopy(tall[tind]);

gind = findall((35.0 .>= temp_fit .>= 0.0) .& (7.0 .>= cond_fit .>= 3.0) .& (50.0 .>= pres_fit .>= 0.0) .& (20.0 .>= chla_fit .>= 0));
dtctdf = dtctd_fit[gind];
tctdf = tctd_fit[gind];
tempf = temp_fit[gind];
condf = cond_fit[gind];
presf = pres_fit[gind];
chlaf = chla_fit[gind];

si = sortperm(tctdf);
tctdf = tctdf[si];
dtctdf = dtctdf[si];
tempf = tempf[si];
condf = condf[si];
presf = presf[si];
chlaf = chlaf[si];
=#

# raw values from the sensor
ttraw = tctd; 
ppraw = presraw[tctdP,2];
zzraw = gsw.gsw_z_from_p.(ppraw*10, llat, 0.0, 0.0); 
ttempraw = tempraw[tctdT,2];
ccondraw = condraw[tctdC,2];
ssaltraw = gsw.gsw_sp_from_c.(ccondraw*10, ttempraw, ppraw*10);
saltAraw= gsw.gsw_sa_from_sp.(ssaltraw, ppraw*10, llon, llat);
ctempraw = gsw.gsw_ct_from_t.(saltAraw, ttempraw, ppraw*10);
rhoraw = gsw.gsw_rho.(saltAraw, ctempraw, ppraw*10);
sigma0raw = gsw.gsw_sigma0.(saltAraw, ctempraw);
spice0raw = gsw.gsw_spiciness0.(saltAraw, ctempraw);
sndspdraw = gsw.gsw_sound_speed.(saltAraw, ctempraw, ppraw*10);

#=
# values fitted to a 1 second time grid
zzf = gsw.gsw_z_from_p.(presf*10, llat, 0.0, 0.0);
saltf = gsw.gsw_sp_from_c.(condf*10, tempf, presf*10);
saltAf = gsw.gsw_sa_from_sp.(saltf, presf*10, llon, llat);
ctempf = gsw.gsw_ct_from_t.(saltAf, tempf, presf*10);
rhof = gsw.gsw_rho.(saltAf, ctempf, presf*10);
sigma0f = gsw.gsw_sigma0.(saltAf, ctempf);
spice0f = gsw.gsw_spiciness0.(saltAf, ctempf);
sndspdf = gsw.gsw_sound_speed.(saltAf, ctempf, presf*10);
=#

#engData = engStruct[];
#sciData = sciStruct[];
ctdData = ctdStruct(mission, glidername, ttraw, ppraw, zzraw, m_lon[:,2], m_lat[:,2], ttempraw, ccondraw, ssaltraw, ctempraw, saltAraw, sigma0raw, spice0raw, sndspdraw, 0, 0);


#x = ttraw;
#y = zzraw;
#z = ttempraw;

x = sci_water_pressure[end-1000:end,1];
y = sci_water_pressure[end-1000:end,2];
z = sci_water_temp[end-1000:end,2];

#x = temptime;
#y = tempz;
#z = tempraw[:,2];

zmin = NaNMath.minimum(z);
zmax = NaNMath.maximum(z); 
fig = Figure(resolution = pres)
ax = Axis(fig[1, 1],
    title = mission * " " * glidername * " Temperature",
    xlabel = "Time",
    ylabel = "Pressure (bar)"
)
Makie.scatter!(x, y, color=z, colormap=:jet, markersize=ms, colorrange=(zmin, zmax))
ax.xticks = (xtick, xticklabel);
Colorbar(fig[1, 2], limits = (zmin, zmax), colormap = :jet, flipaxis = false)
fig
save(figoutdir * mission * "_" * glidername * "_temp.png", fig)

=#
=#