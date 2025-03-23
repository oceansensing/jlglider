using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations, YAML

using Glider
import Glider.slocumType: engStruct, ctdStruct, sciStruct

include("/Users/gong/GitHub/jlglider/common/C2PO.jl")
import .C2PO: pyrow2jlcol, intersectalajulia2, intersectalajulia4, unix2yearday, yearday2unix, datetime2yearday, yearday2datetime

include("slocumFunc.jl")
using .slocumFunc: glider_var_load, glider_presfunc

missionYAMLdirpath = "/Users/gong/GitHub/jlglider/slocum/slocum_yaml_PASSENGERS/2024/";

if missionYAMLdirpath[end-4:end] != ".yaml"
    if missionYAMLdirpath[end] != '/'
        missionYAMLdirpath *= "/";
    end
    missionYAMLpath = Glob.glob("*.yaml", missionYAMLdirpath);
else
    missionYAMLpath = [missionYAMLdirpath];
end

i = 2

dbdreader = pyimport("dbdreader");
gsw = GibbsSeaWater;

tbound = 3600*24*365.0 * 0.2; 
latmin, latmax = 38, 40;
lonmin, lonmax = -65.2, -61.3;    

# load mission YAML file
missionYAML = YAML.load(IOBuffer(read(missionYAMLpath[i], String)));
glidertype = missionYAML["gliderType"];
glidername = missionYAML["gliderName"];
gliderSN = missionYAML["gliderSN"];
missionID = missionYAML["missionID"];
project = missionYAML["project"];
dataroot = missionYAML["dataroot"];
deploydate = string(missionYAML["deploydate"]);
datamode = missionYAML["dataflag"];
suffix = missionYAML["suffix"];

datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/"; 
dataENGdir = datadir * "eng/";
dataSCIdir = datadir * "sci/";
cacdir = datadir * "cache/";

if datamode == "realtime"
    dataGliderEng = dbdreader.MultiDBD(pattern = dataENGdir * "*.[sS][bB][dD]", cacheDir = cacdir, skip_initial_line = true);
    dataGliderSci = dbdreader.MultiDBD(pattern = dataSCIdir * "*.[tT][bB][dD]", cacheDir = cacdir, skip_initial_line = true);
else
    dataGliderEng = dbdreader.MultiDBD(pattern = dataENGdir * "*.[dD][bB][dD]", cacheDir = cacdir, skip_initial_line = true);
    dataGliderSci = dbdreader.MultiDBD(pattern = dataSCIdir * "*.[eE][bB][dD]", cacheDir = cacdir, skip_initial_line = true);
end

engvars = dataGliderEng.parameterNames["eng"];
scivars = dataGliderSci.parameterNames["sci"];

# unix time of the mission start time
unixt0 = Dates.datetime2unix(Dates.DateTime(deploydate[1:4] * "-" * deploydate[5:6] * "-" * deploydate[7:8]));
trange = [unixt0, unixt0 + tbound]; # mission start time + 6 months

# load engineering and CTD data from raw glider DBD and EBD files
m_present_time = dataGliderEng.get("m_present_time")[1];
m_present_time_ind = findall(unixt0 .< m_present_time .< (unixt0 + tbound));
#m_present_time_ind = findall(trange[1] .< m_present_time .< trange[end]);
m_present_time = m_present_time[m_present_time_ind];

sci_m_present_time = dataGliderEng.get("sci_m_present_time")[1];
sci_m_present_time_ind = findall(unixt0 .< sci_m_present_time .< (unixt0 + tbound));
#sci_m_present_time_ind = findall(trange[1] .< sci_m_present_time .< trange[end]);
sci_m_present_time = sci_m_present_time[sci_m_present_time_ind];

m_gps_lat = dataGliderEng.get("m_gps_lat", return_nans=false);
m_gps_lon = dataGliderEng.get("m_gps_lon", return_nans=false);

midlat = median(m_gps_lat[2]);
midlon = median(m_gps_lon[2]);
latmin, latmax = midlat - 10, midlat + 10;
lonmin, lonmax = midlon - 10, midlon + 10;    

latind = findall((latmin .<= m_gps_lat[2] .<= latmax) .& (1e9 .< m_gps_lat[1] .< 1e10));
m_gps_lat_t = m_gps_lat[1][latind]
m_gps_lat_v = m_gps_lat[2][latind]
m_gps_lat = (m_gps_lat_t, m_gps_lat_v);

lonind = findall((lonmin .<= m_gps_lon[2] .<= lonmax) .& (1e9 .< m_gps_lon[1] .< 1e10));
m_gps_lon_t = m_gps_lon[1][lonind]
m_gps_lon_v = m_gps_lon[2][lonind]
m_gps_lon = (m_gps_lon_t, m_gps_lon_v);

sci_water_pressure = dataGliderSci.get("sci_water_pressure", return_nans=true);
sci_water_temp = dataGliderSci.get("sci_water_temp", return_nans=true);
sci_water_cond = dataGliderSci.get("sci_water_cond", return_nans=true);

mlon = NaNMath.mean(m_gps_lon[2]);
mlat = NaNMath.mean(m_gps_lat[2]);

presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, trange);
lonfunc, lontime, lonraw, lonpres, lonz = glider_var_load(m_gps_lon, trange, [-80.0 10.0], presfunc, mlat);
latfunc, lattime, latraw, latpres, latz = glider_var_load(m_gps_lat, trange, [20.0 80.0], presfunc, mlat);  
tempfunc, temptime, tempraw, temppres, tempz = glider_var_load(sci_water_temp, trange, [-2.0 40.0], presfunc, mlat)
condfunc, condtime, condraw, condpres, condz = glider_var_load(sci_water_cond, trange, [0.1 100.0], presfunc, mlat)

#presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, trange);
#lonfunc, lontime, lonraw, lonpres, lonz = glider_var_load(m_gps_lon, trange, [-80.0 -50.0], sci_water_pressure, mlat);
#latfunc, lattime, latraw, latpres, latz = glider_var_load(m_gps_lat, trange, [20.0 60.0], sci_water_pressure, mlat);  
#tempfunc, temptime, tempraw, temppres, tempz = glider_var_load(sci_water_temp, trange, [0.1 40.0], sci_water_pressure, mlat)
#condfunc, condtime, condraw, condpres, condz = glider_var_load(sci_water_cond, trange, [0.01 100.0], sci_water_pressure, mlat)


# find common glider values
tctd = unique(intersect(prestime, temptime, condtime));
tctdT = intersectalajulia2(tctd, temptime)[3];
tctdC = intersectalajulia2(tctd, condtime)[3];
tctdP = intersectalajulia2(tctd, prestime)[3];

lonf = lonfunc(tctd);
latf = latfunc(tctd);

pres = presraw[tctdP];
temp = tempraw[tctdT];
cond = condraw[tctdC];

# apply Gibbs SeaWater functions to calculate derived values
z = gsw.gsw_z_from_p.(pres*10, mlat, 0.0, 0.0); 
salt = gsw.gsw_sp_from_c.(cond*10, temp, pres*10);
saltA = gsw.gsw_sa_from_sp.(salt, pres*10, mlon, mlat);
ctemp = gsw.gsw_ct_from_t.(saltA, temp, pres*10);
rho = gsw.gsw_rho.(saltA, ctemp, pres*10);
sigma0 = gsw.gsw_sigma0.(saltA, ctemp);
spice0 = gsw.gsw_spiciness0.(saltA, ctemp);
sndspd = gsw.gsw_sound_speed.(saltA, ctemp, pres*10);

# save glider CTD data to a ctdStruct object
ctdData = Glider.slocumType.ctdStruct(glidertype, gliderSN, glidername, missionID, project, tctd, pres, z, lonf, latf, temp, cond, salt, ctemp, saltA, sigma0, spice0, sndspd, 0, 0);