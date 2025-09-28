# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
# gong@vims.edu 2024-11-11: add slocumYAMLload function and load_glider_ctd using YAML metadata

module slocumLoad

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations, YAML

using Glider
import Glider.slocumType: engStruct, ctdStruct, sciStruct

include("/Users/gong/GitHub/jlglider/common/C2PO.jl")
import .C2PO: pyrow2jlcol, intersectalajulia2, intersectalajulia4, unix2yearday, yearday2unix, datetime2yearday, yearday2datetime

include("slocumFunc.jl")
using .slocumFunc: glider_var_load, glider_presfunc

function glider_ctd_qc(t::Array{Float64}, pres::Array{Float64}, temp::Array{Float64}, cond::Array{Float64}, trange::Array{Float64})
    gind = findall((trange[1] .<= t .<= trange[end]) .& (0.0 .<= pres .<= 105.0) .& (-2.0 .<= temp .<= 40.0) .& (2.0 .<= cond .<= 7.0)); 
    return t[gind], pres[gind], temp[gind], cond[gind];
end

function load_glider_ctd(datadir, cacdir, trange, lonrange, latrange, datamode, glidertype, gliderSN, glidername, missionID, project, loadmode)
    dbdreader = pyimport("dbdreader");
    gsw = GibbsSeaWater;
    
    engdatadir = datadir * "gli/";
    scidatadir = datadir * "sci/";

    if loadmode == "uppercase"
        if datamode == "realtime"
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = engdatadir * "*.[sStT][bB][dD]", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true);
        else
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = scidatadir * "*.[dDeE][bB][dD]", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true);
        end
    else
        if datamode == "realtime"
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = engdatadir * "*.[sStT][bB][dD]", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true);
        else
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = scidatadir * "*.[dDeE][bB][dD]", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true);
        end
    end

    engvars = dataGlider.parameterNames["eng"];
    scivars = dataGlider.parameterNames["sci"];

    # load engineering data from raw glider DBD files
    m_present_time = dataGlider.get("m_present_time");
    m_lat = dataGlider.get("m_lat"); 
    m_lon = dataGlider.get("m_lon"); 
    m_gps_lat = dataGlider.get("m_gps_lat");
    m_gps_lon = dataGlider.get("m_gps_lon");
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
    sci_m_present_time = dataGlider.get("sci_m_present_time");
    sci_water_temp = dataGlider.get("sci_water_temp");
    sci_water_cond = dataGlider.get("sci_water_cond");
    sci_water_pressure = dataGlider.get("sci_water_pressure");
    #sci_flbbcd_chlor_units = dataGlider.get("sci_flbbcd_chlor_units");
    #sci_flbbcd_cdom_units = dataGlider.get("sci_flbbcd_cdom_units");
    #sci_flbbcd_bb_units = dataGlider.get("sci_flbbcd_bb_units");
    #sci_bsipar_par = dataGlider.get("sci_bsipar_par");

    # set lon & lat limits and calculate mean lon and lat values.
    #glonind = findall(-65 .< m_lon[2,:] .< -55); 
    #glatind = findall(30 .< m_lat[2,:] .< 50);
    #lonrange = [-65.0 -55.0];
    #latrange = [30.0 50.0]; 
    glonind = findall(lonrange[1] .< m_lon[2] .< lonrange[2]);
    glatind = findall(latrange[1] .< m_lat[2] .< latrange[2]);  
    mlat = Statistics.mean(m_lat[2][glatind]);
    mlon = Statistics.mean(m_lon[2][glonind]);    #llon = -73.4;
    #llat = 38.0;

    # calculate derived values from CTD data
    presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, trange);
    tempfunc, temptime, temppres, tempraw, tempz = glider_var_load(sci_water_temp, trange, [0.1 40.0], sci_water_pressure, mlat)
    condfunc, condtime, condpres, condraw, condz = glider_var_load(sci_water_cond, trange, [0.01 100.0], sci_water_pressure, mlat)

    lonfunc, lontime, lonpres, lonraw, lonz = glider_var_load(m_gps_lon, trange, lonrange, sci_water_pressure, mlat);
    latfunc, lattime, latpres, latraw, latz = glider_var_load(m_gps_lat, trange, latrange, sci_water_pressure, mlat); 

    # find common glider values
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
    salt = gsw.gsw_sp_from_c.(cond*10, temp, pres*10);
    saltA= gsw.gsw_sa_from_sp.(salt, pres*10, mlon, mlat);
    ctemp = gsw.gsw_ct_from_t.(saltA, temp, pres*10);
    rho = gsw.gsw_rho.(saltA, ctemp, pres*10);
    sigma0 = gsw.gsw_sigma0.(saltA, ctemp);
    spice0 = gsw.gsw_spiciness0.(saltA, ctemp);
    sndspd = gsw.gsw_sound_speed.(saltA, ctemp, pres*10);

    ctdData = Glider.slocumType.ctdStruct(glidertype, gliderSN, glidername, missionID, project, tctd, pres, z, lonf, latf, temp, cond, salt, ctemp, saltA, sigma0, spice0, sndspd, 0, 0);
    return ctdData
end

# new function to load glider CTD data using YAML metadata file (2024-09-05)
function load_glider_ctd(missionYAMLpath::String)
    dbdreader = pyimport("dbdreader");
    gsw = GibbsSeaWater;

    tbound = 3600*24*365.0 * 0.2; 

    # load mission YAML file
    missionYAML = YAML.load(IOBuffer(read(missionYAMLpath, String)));
    glidertype = missionYAML["gliderType"];
    glidername = missionYAML["gliderName"];
    gliderSN = missionYAML["gliderSN"];
    missionID = missionYAML["missionID"];
    project = missionYAML["project"];
    dataroot = missionYAML["dataroot"];
    deploydate = string(missionYAML["deploydate"]);
    datamode = missionYAML["dataflag"];
    suffix = missionYAML["suffix"];

    #datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/" * "data/"; 
    datadir = dataroot * glidername * "-" * deploydate * "-" * project * "/" * "data/"; 
    dataENGdir = datadir * "gli/";
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
    lonmin, lonmax = midlon - 20, midlon + 20;    

    latind = findall((latmin .<= m_gps_lat[2] .<= latmax) .& (1e9 .< m_gps_lat[1] .< 1e10));
    m_gps_lat_t = m_gps_lat[1][latind];
    m_gps_lat_v = m_gps_lat[2][latind];
    m_gps_lat = (m_gps_lat_t, m_gps_lat_v);

    lonind = findall((lonmin .<= m_gps_lon[2] .<= lonmax) .& (1e9 .< m_gps_lon[1] .< 1e10));
    m_gps_lon_t = m_gps_lon[1][lonind];
    m_gps_lon_v = m_gps_lon[2][lonind];
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
    #[Tu, Rsubrho, p_mid] = gsw.gsw_Turner_Rsubrho(saltA,ctemp,pres*10);

    # save glider CTD data to a ctdStruct object
    ctdData = Glider.slocumType.ctdStruct(glidertype, gliderSN, glidername, missionID, project, tctd, pres, z, lonf, latf, temp, cond, salt, ctemp, saltA, sigma0, spice0, sndspd, 0, 0);
    return ctdData
end

function load_glider_ctd(datadir, cacdir, trange, datamode, mission, glidername)
    dbdreader = pyimport("dbdreader");
    gsw = GibbsSeaWater;

    gliderSN = -1;
    
    if (@isdefined datadir) == false
        # setup directories
        rootdir = "/Users/gong/oceansensing Dropbox/C2PO/MARACOOS/";
        fromgliderdir = rootdir * "electa-20230320-maracoos/from-glider/"; 
        datadirpath = Glob.glob("electa-from-glider*.zip", fromgliderdir)[1];
        datadir = fromgliderdir * "electa-from-glider-20230404T113550/";
        cacdir = rootdir * "electa-20230320-maracoos/from-glider/cache/";
    end

    if (@isdefined trange) == false
        # specify valid data time period
        t0 = DateTime("2023-03-21");
        tN = DateTime("2023-04-21");
        trange = datetime2unix.([t0; tN]);
    end

    # setup glider data loading using dbdreader
   if datamode == "realtime"
        #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
        dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[sStT]*", cacheDir = cacdir, complement_files_only = true, skip_initial_line = true);
    else
        #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files = true, cacheDir = cacdir);
        dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[dDeE]*", cacheDir = cacdir, complement_files_only = true, skip_initial_line = true);
    end
    engvars = dataGlider.parameterNames["eng"];
    scivars = dataGlider.parameterNames["sci"];

    t, sci_m_present_time, lon, lat, m_pressure, sci_water_pressure, sci_water_temp, sci_water_cond = dataGlider.get_sync("sci_m_present_time", "m_lon", "m_lat", "m_pressure", "sci_water_pressure", "sci_water_temp", "sci_water_cond");
    #t2, sci_m_present_time2, sci_flbbcd_chlor_units, sci_flbbcd_cdom_units, sci_flbbcd_bb_units, sci_bsipar_par = dataGlider.get_sync("sci_m_present_time", "sci_flbbcd_chlor_units", "sci_flbbcd_cdom_units", "sci_flbbcd_bb_units", "sci_bsipar_par");
    tis = sortperm(sci_m_present_time);

    mlon = NaNMath.mean(lon);
    mlat = NaNMath.mean(lat);

    tctd, pres, temp, cond = glider_ctd_qc(sci_m_present_time[tis], sci_water_pressure[tis], sci_water_temp[tis], sci_water_cond[tis], trange);
    #tctd2, pres2, temp2, temp2ind = glider_var_load(sci_m_present_time, sci_water_pressure, sci_water_temp, trange, [0.1 40.0]);

    z = gsw.gsw_z_from_p.(pres*10, mlat, 0.0, 0.0); 
    salt = gsw.gsw_sp_from_c.(cond*10, temp, pres*10);
    saltA = gsw.gsw_sa_from_sp.(salt, pres*10, mlon, mlat);
    ctemp = gsw.gsw_ct_from_t.(saltA, temp, pres*10);
    rho = gsw.gsw_rho.(saltA, ctemp, pres*10);
    sigma0 = gsw.gsw_sigma0.(saltA, ctemp);
    spice0 = gsw.gsw_spiciness0.(saltA, ctemp);
    sndspd = gsw.gsw_sound_speed.(saltA, ctemp, pres*10);

    
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
    sci_flbbcd_chlor_units = dataGlider.get("sci_flbbcd_chlor_units");
    sci_flbbcd_cdom_units = dataGlider.get("sci_flbbcd_cdom_units");
    sci_flbbcd_bb_units = dataGlider.get("sci_flbbcd_bb_units");
    sci_bsipar_par = dataGlider.get("sci_bsipar_par");
    =#

    #=
    # calculate derived values from CTD data
    llat = Statistics.mean(m_lat[2]);
    llon = Statistics.mean(m_lon[2]);
    #llon = -73.4;
    #llat = 38.0;

    presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, trange);
    tempraw, temptime, temppres, tempz = glider_var_load(sci_water_temp, trange, [0.1 40.0], sci_water_pressure, llat)
    condraw, condtime, condpres, condz = glider_var_load(sci_water_cond, trange, [0.01 100.0], sci_water_pressure, llat)

    # find common glider values
    tctd = unique(intersect(prestime, temptime, condtime));
    tctdT = intersectalajulia2(tctd, temptime)[3];
    tctdC = intersectalajulia2(tctd, condtime)[3];
    tctdP = intersectalajulia2(tctd, prestime)[3];

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
    ppraw = presraw[tctdP];
    zzraw = gsw.gsw_z_from_p.(ppraw*10, llat, 0.0, 0.0); 
    ttempraw = tempraw[tctdT];
    ccondraw = condraw[tctdC];
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
    =#

    #engData = engStruct[];
    #sciData = sciStruct[];
    ctdData = Glider.slocumType.ctdStruct(mission, gliderSN, glidername, tctd, pres, z, lon, lat, temp, cond, salt, ctemp, saltA, sigma0, spice0, sndspd, 0, 0);
    return ctdData
end

function load_glider_sci(datadir, cacdir, trange, datamode, mission, glidername, loadmode)
    dbdreader = pyimport("dbdreader");
    gsw = GibbsSeaWater;

    # setup glider data loading using dbdreader
    if datamode == "realtime"
        #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
        dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", cacheDir = cacdir, complement_files_only = true, skip_initial_line = true);
    else
        dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", cacheDir = cacdir, complement_files_only = true, skip_initial_line = true);
    end
    engvars = dataGlider.parameterNames["eng"];
    scivars = dataGlider.parameterNames["sci"];

    # load engineering data from raw glider DBD files
    m_present_time = dataGlider.get("m_present_time");
    m_lat = dataGlider.get("m_lat"); 
    m_lon = dataGlider.get("m_lon"); 

    # load CTD data from raw glider EBD files
    sci_m_present_time = dataGlider.get("sci_m_present_time");
    sci_water_temp = dataGlider.get("sci_water_temp");
    sci_water_cond = dataGlider.get("sci_water_cond");
    sci_water_pressure = dataGlider.get("sci_water_pressure");
    sci_flbbcd_chlor_units = dataGlider.get("sci_flbbcd_chlor_units");
    sci_flbbcd_cdom_units = dataGlider.get("sci_flbbcd_cdom_units");
    sci_flbbcd_bb_units = dataGlider.get("sci_flbbcd_bb_units");
    sci_bsipar_par = dataGlider.get("sci_bsipar_par");

    # calculate derived values from CTD data
    mlat = Statistics.mean(m_lat[2]);
    mlon = Statistics.mean(m_lon[2]);

    presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, trange);
    tempfunc, temptime, temppres, tempraw, tempz = glider_var_load(sci_water_temp, trange, [0.1 40.0], sci_water_pressure, mlat)
    condfunc, condtime, condpres, condraw, condz = glider_var_load(sci_water_cond, trange, [0.01 100.0], sci_water_pressure, mlat)

    # find common glider values
    tctd = unique(intersect(prestime, temptime, condtime));
    tind = findall(trange[1] .<= tctd .<= trange[end]);

    #t, sci_m_present_time, lon, lat, sci_water_pressure, sci_flbbcd_chlor_units, sci_flbbcd_cdom_units, sci_flbbcd_bb_units, sci_bsipar_par = dataGlider.get_sync("sci_m_present_time", "m_lon", "m_lat", "sci_water_pressure", "sci_flbbcd_chlor_units", "sci_flbbcd_cdom_units", "sci_flbbcd_bb_units", "sci_bsipar_par");
    #tis = sortperm(sci_m_present_time);

    #mlon = NaNMath.mean(lon);
    #mlat = NaNMath.mean(lat);

    #=
    m_lat = dataGlider.get("m_lat"); 
    m_lon = dataGlider.get("m_lon"); 

    sci_m_present_time = dataGlider.get("sci_m_present_time");
    sci_water_pressure = dataGlider.get("sci_water_pressure");
    sci_flbbcd_chlor_units = dataGlider.get("sci_flbbcd_chlor_units");
    sci_flbbcd_cdom_units = dataGlider.get("sci_flbbcd_cdom_units");
    sci_flbbcd_bb_units = dataGlider.get("sci_flbbcd_bb_units");
    sci_bsipar_par = dataGlider.get("sci_bsipar_par");

    llat = Statistics.mean(m_lat[2]);
    llon = Statistics.mean(m_lon[2]);
    =#

    #sci_m_present_time, sci_water_pressure, sci_water_temp, trange, [0.1 40.0]

    if isempty(sci_flbbcd_chlor_units) != true
        chlatime, chlapres, chlaraw, chlaind = glider_var_load(sci_m_present_time[tis], sci_water_pressure[tis], sci_flbbcd_chlor_units, trange, [-0.1 3.0])
        chlaz = gsw.gsw_z_from_p.(chlapres*10, mlat, 0.0, 0.0);
        chlaData = Glider.slocumType.sciStruct(mission, glidername, chlatime, chlapres, chlaz, lon, lat, chlaraw);
    else
        chlaData = [];
    end

    if isempty(sci_flbbcd_cdom_units) != true
        cdomtime, cdompres, cdomraw, cdomind = glider_var_load(sci_m_present_time, sci_water_pressure, sci_flbbcd_cdom_units, trange, [-5.0 5.0])
        cdomz = gsw.gsw_z_from_p.(cdompres*10, mlat, 0.0, 0.0); 
        cdomData = Glider.slocumType.sciStruct(mission, glidername, cdomtime, cdompres, cdomz, lon, lat, cdomraw);
    else
        cdomData = [];
    end

    if isempty(sci_flbbcd_bb_units) != true
        bb700time, bb700pres, bb700raw, bb700ind = glider_var_load(sci_m_present_time, sci_water_pressure, sci_flbbcd_bb_units, trange, [0.0 0.008])
        bb700z = gsw.gsw_z_from_p.(bb700pres*10, mlat, 0.0, 0.0); 
        bb700Data = Glider.slocumType.sciStruct(mission, glidername, bb700time, bb700pres, bb700z, lon, lat, bb700raw);
    else
        bb700Data = [];
    end

    if isempty(sci_bsipar_par) != true
        bpartime, bparpres, bparraw, bparind = glider_var_load(sci_m_present_time, sci_water_pressure, sci_bsipar_par, trange, [0.0 6000.0])
        bparz = gsw.gsw_z_from_p.(bparpres*10, mlat, 0.0, 0.0); 
        bparData = Glider.slocumType.sciStruct(mission, glidername, bpartime, bparpres, bparz, lon, lat, bparraw);
    else
        bparData = [];
    end

    return chlaData, cdomData, bb700Data, bparData;
end

# this function load the glider CTD data from the gliderData directory using metadata from YAML file for each mission
function slocumYAMLload(missionYAMLdirpath::String)
    if (@isdefined missionYAMLdirpath) == false
        missionYAMLdirpath = "/Users/gong/GitHub/jlglider/slocum/mission_yaml_GS/";
    end

    if missionYAMLdirpath[end-4:end] != ".yaml"
        if missionYAMLdirpath[end] != '/'
         missionYAMLdirpath *= "/";
        end
        missionYAMLpath = Glob.glob("*.yaml", missionYAMLdirpath);
    else
        missionYAMLpath = [missionYAMLdirpath];
    end

    gliderCTDarray = ctdStruct[];

    for i = 1:length(missionYAMLpath)
        display(missionYAMLpath[i])
        push!(gliderCTDarray, load_glider_ctd(missionYAMLpath[i]));
    end
    return gliderCTDarray;
end

end