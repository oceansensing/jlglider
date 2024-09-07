using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations, YAML

include("slocumType.jl")
include("slocumFunc.jl")

#import .slocumType: ctdStruct, sciStruct
using .slocumType: ctdStruct, sciStruct
using .slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_presfunc, yearday2datetime, datetime2yearday

missionYAMLpath = "/Users/gong/GitHub/jlglider/slocum/mission_yaml/electa-20230321-maracoos.yaml"

dbdreader = pyimport("dbdreader");
gsw = GibbsSeaWater;

# load mission YAML file
missionYAML = YAML.load(IOBuffer(read(missionYAMLpath, String)));
glidertype = missionYAML["gliderType"];
glidername = missionYAML["gliderName"];
gliderSN = missionYAML["gliderSN"];
project = missionYAML["project"];
dataroot = missionYAML["dataroot"];
deploydate = string(missionYAML["deploydate"]);
datamode = missionYAML["dataflag"];
suffix = missionYAML["suffix"];

datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
cacdir = datadir * "cache/";

if datamode == "realtime"
    #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
    dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[sStT][bB][dD]", cacheDir = cacdir, complemented_files_only = true, skip_initial_line = true);
else
    #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files = true, cacheDir = cacdir);
    dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[dDeE][bB][dD]", cacheDir = cacdir, complemented_files_only = true, skip_initial_line = true);
end

engvars = dataGlider.parameterNames["eng"];
scivars = dataGlider.parameterNames["sci"];

swt = dataGlider.get("sci_water_temp", decimalLatLon=true, discardBadLatLon=true, return_nans=true, include_source=false, max_values_to_read=-1)[1]
gswt = findall(1.0 .< swt .< 40)

tbound = 3600*24*365.0;

# load engineering and CTD data from raw glider DBD and EBD files
m_present_time = dataGlider.get("m_present_time")[1];
m_present_time_ind = findall(NaNMath.median(m_present_time) - tbound .< m_present_time .< NaNMath.median(m_present_time) + tbound);
m_present_time = m_present_time[m_present_time_ind];

sci_m_present_time = dataGlider.get("sci_m_present_time")[1];
sci_m_present_time_ind = findall(NaNMath.median(sci_m_present_time) - tbound .< sci_m_present_time .< NaNMath.median(sci_m_present_time) + tbound);
sci_m_present_time = sci_m_present_time[sci_m_present_time_ind];

m_gps_lat = dataGlider.get("m_gps_lat", return_nans=false);
#m_lat_t = m_lat[1];
#m_lat = m_lat[2];
#m_lat_ind = findall((median(m_lat_t) - tbound .< m_lat_t .< median(m_lat_t) + tbound) .& (0.1 .< m_lat .< 90.0));
#m_lat_t = m_lat_t[m_lat_ind];
#m_lat = m_lat[m_lat_ind];

m_gps_lon = dataGlider.get("m_gps_lon", return_nans=false); 
#m_lon_t = m_lon[1];
#m_lon = m_lon[2];
#m_lon_ind = findall((median(m_lon_t) - tbound .< m_lon_t .< median(m_lon_t) + tbound) .& (-180.0 .< m_lon .< 180.0));
#m_lon_t = m_lon_t[m_lon_ind];
#m_lon = m_lon[m_lon_ind];

sci_water_pressure = dataGlider.get("sci_water_pressure", return_nans=true);
#sci_water_pressure_t = sci_water_pressure[1];
#sci_water_pressure = sci_water_pressure[2];
#sci_water_pressure_ind = findall((median(sci_water_pressure_t) - tbound .< sci_water_pressure_t .< median(sci_water_pressure_t) + tbound) .& (0.0 .< sci_water_pressure .< 105.0));
#sci_water_pressure_t = sci_water_pressure_t[sci_water_pressure_ind];
#sci_water_pressure = sci_water_pressure[sci_water_pressure_ind];

sci_water_temp = dataGlider.get("sci_water_temp", return_nans=true);
#sci_water_temp_t = sci_water_temp[1];
#sci_water_temp = sci_water_temp[2];
#sci_water_temp_ind = findall((median(sci_water_temp_t) - tbound .< sci_water_temp_t .< median(sci_water_temp_t) + tbound) .& (0.1 .< sci_water_temp .< 40.0));
#sci_water_temp_t = sci_water_temp_t[sci_water_temp_ind];
#sci_water_temp = sci_water_temp[sci_water_temp_ind];

sci_water_cond = dataGlider.get("sci_water_cond", return_nans=true);
#sci_water_cond_t = sci_water_cond[1];
#sci_water_cond = sci_water_cond[2];
#sci_water_cond_ind = findall((median(sci_water_cond_t) - tbound .< sci_water_cond_t .< median(sci_water_cond_t) + tbound) .& (2.0 .< sci_water_cond .< 7.0));
#sci_water_cond_t = sci_water_cond_t[sci_water_cond_ind];
#sci_water_cond = sci_water_cond[sci_water_cond_ind];


#t, sci_m_present_time, lon, lat, m_pressure, sci_water_pressure, sci_water_temp, sci_water_cond = dataGlider.get_sync("sci_m_present_time", "m_lon", "m_lat", "m_pressure", "sci_water_pressure", "sci_water_temp", "sci_water_cond");
#t2, sci_m_present_time2, sci_flbbcd_chlor_units, sci_flbbcd_cdom_units, sci_flbbcd_bb_units, sci_bsipar_par = dataGlider.get_sync("sci_m_present_time", "sci_flbbcd_chlor_units", "sci_flbbcd_cdom_units", "sci_flbbcd_bb_units", "sci_bsipar_par");

mlon = NaNMath.mean(m_gps_lon[2]);
mlat = NaNMath.mean(m_gps_lat[2]);

presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, tbound);

lonfunc, lontime, lonpres, lonraw, lonz = glider_var_load(m_gps_lon, tbound, [-80.0 -50.0], presfunc, mlat);
latfunc, lattime, latpres, latraw, latz = glider_var_load(m_gps_lat, tbound, [20.0 60.0], presfunc, mlat);

tempfunc, temptime, temppres, tempraw, tempz = glider_var_load(sci_water_temp, tbound, [0.1 40.0], presfunc, mlat)
condfunc, condtime, condpres, condraw, condz = glider_var_load(sci_water_cond, tbound, [0.01 100.0], presfunc, mlat)

glidervar = sci_water_temp;
varlim = [0.1 40.0];

varind = findall(((NaNMath.median(glidervar[1]) - tbound) .<= glidervar[1] .<= (NaNMath.median(glidervar[1]) + tbound)) .& (varlim[1] .<= glidervar[2] .<= varlim[end])); 
vartime = glidervar[1][varind];
varraw = glidervar[2][varind];
sortedvarind = sortperm(vartime);
varrawval = varraw[sortedvarind];
vartime = vartime[sortedvarind];
varfunc = linear_interpolation(vartime, varrawval, extrapolation_bc=Line()); 
vardtime = unix2datetime.(vartime);
varpres = presfunc(vartime);


tctd = unique(intersect(prestime, temptime, condtime));
tctdT = intersectalajulia2(tctd, temptime)[3];
tctdC = intersectalajulia2(tctd, condtime)[3];
tctdP = intersectalajulia2(tctd, prestime)[3];

lonf = lonfunc(tctd);
latf = latfunc(tctd);

pres = presraw[tctdP];
z = gsw.gsw_z_from_p.(pres*10, mlat, 0.0, 0.0); 
temp = tempraw[tctdT];
cond = condraw[tctdC];

#=
# this step set the glider mission's time limits to +/- 365 days of the median time value, should apply to most glider missions except for extremely long ones.
tmedian = median(t);
trange = [tmedian - 3600*24*365; tmedian + 3600*24*365];
#tind = findall(trange[1] .< t .< trange[2]);
gind = findall((trange[1] .<= sci_m_present_time .<= trange[end]) .& (0.0 .<= sci_water_pressure .<= 105.0) .& (-2.0 .<= sci_water_temp .<= 40.0) .& (2.0 .<= sci_water_cond .<= 7.0)); 

t = t[gind];
sci_m_present_time = sci_m_present_time[gind];
lon = lon[gind];
lat = lat[gind];
m_pressure = m_pressure[gind];
sci_water_pressure = sci_water_pressure[gind];
sci_water_temp = sci_water_temp[gind];
sci_water_cond = sci_water_cond[gind];

tis = sortperm(sci_m_present_time);

tctd, pres, temp, cond = sci_m_present_time[tis], sci_water_pressure[tis], sci_water_temp[tis], sci_water_cond[tis];

#tctd, pres, temp, cond = glider_ctd_qc(sci_m_present_time[tis], sci_water_pressure[tis], sci_water_temp[tis], sci_water_cond[tis], trange);
#tctd2, pres2, temp2, temp2ind = glider_var_load(sci_m_present_time, sci_water_pressure, sci_water_temp, trange, [0.1 40.0]);
=#


#z = gsw.gsw_z_from_p.(pres*10, mlat, 0.0, 0.0); 
salt = gsw.gsw_sp_from_c.(cond*10, temp, pres*10);
saltA = gsw.gsw_sa_from_sp.(salt, pres*10, mlon, mlat);
ctemp = gsw.gsw_ct_from_t.(saltA, temp, pres*10);
rho = gsw.gsw_rho.(saltA, ctemp, pres*10);
sigma0 = gsw.gsw_sigma0.(saltA, ctemp);
spice0 = gsw.gsw_spiciness0.(saltA, ctemp);
sndspd = gsw.gsw_sound_speed.(saltA, ctemp, pres*10);

ctdData = ctdStruct(project, glidername, tctd, pres, z, lonf, latf, temp, cond, salt, ctemp, saltA, sigma0, spice0, sndspd, 0, 0);