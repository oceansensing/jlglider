# slocumLoadTest

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations

import slocumType: ctdStruct, sciStruct
import slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_ctd_load, glider_presfunc

#dbdreaderdir = "/Users/gong/GitHub/dbdreader/";
#pushfirst!(pyimport("sys")."path", dbdreaderdir);

#function load_glider_ctd(datadir, cacdir, trange, datamode, mission, glidername, loadmode)
    dbdreader = pyimport("dbdreader");
    gsw = GibbsSeaWater;

    if loadmode == "uppercase"
        if datamode == "realtime"
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[ST]BD", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true, decimalLatLon = true, return_nans = true);
        else
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[DE]BD", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true, decimalLatLon = true, return_nans = true);
        end
    else
        if datamode == "realtime"
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true, decimalLatLon = true, return_nans = true);
        else
            #dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files = true, cacheDir = cacdir);
            dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", cacheDir = cacdir, complemented_files_only = false, skip_initial_line = true, decimalLatLon = true, return_nans = true);
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

    # calculate derived values from CTD data
    glonind = findall(-65 .< m_lon[2,:] .< -55); 
    glatind = findall(30 .< m_lat[2,:] .< 50);
    mlat = Statistics.mean(m_lat[2,glatind]);
    mlon = Statistics.mean(m_lon[2,glonind]);
    #llon = -73.4;
    #llat = 38.0;

    #t, C, T, P, mlon, mlat, mpresenttime = dataGlider.get_CTD_sync("m_lon", "m_lat", "m_present_time")

    presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, trange);
    tempfunc, temptime, temppres, tempraw, tempz = glider_var_load(sci_water_temp, trange, [0.1 40.0], sci_water_pressure, mlat)
    condfunc, condtime, condpres, condraw, condz = glider_var_load(sci_water_cond, trange, [0.01 100.0], sci_water_pressure, mlat)

    lonfunc, lontime, lonpres, lonraw, lonz = glider_var_load(m_gps_lon, trange, [-65.0 -55.0], sci_water_pressure, mlat);
    latfunc, lattime, latpres, latraw, latz = glider_var_load(m_gps_lat, trange, [30.0 50.0], sci_water_pressure, mlat); 

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

    ctdData = ctdStruct(mission, glidername, tctd, pres, z, lonf, latf, temp, cond, salt, ctemp, saltA, sigma0, spice0, sndspd, 0, 0);
 