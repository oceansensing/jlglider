# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#

module load_slocum_glider

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations

import slocumType: ctdStruct
import slocumFunc: pyrow2jlcol, intersectalajulia2, glider_var_load, glider_presfunc

#=
# define function for converting an array from python row major to julia column major
function pyrow2jlcol(invar::Matrix{Float64})
    return reverse(rotr90(invar), dims = 2);
end

# https://discourse.julialang.org/t/indices-of-intersection-of-two-arrays/23043/20
function intersectalajulia2(a,b)
    ia = findall(in(b), a)
    ib = findall(in(view(a,ia)), b)
    return unique(view(a,ia)), ia, ib[indexin(view(a,ia), view(b,ib))]
end

function intersectalajulia4(a,b)
    ab=intersect(a,b)
    ia = [findall(==(e), a) for e in ab]
    ib = [findall(==(e), b) for e in ab]
    return hcat(ab, ia,ib)
end
=#

function load_glider_ctd(datadir, cacdir, trange, datamode, mission, glidername)
    dbdreader = pyimport("dbdreader");
    gsw = GibbsSeaWater;
    
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
        dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[st]bd", complement_files = true, cacheDir = cacdir);
    else
        dataGlider = dbdreader.MultiDBD(pattern = datadir * "*.[de]bd", complement_files = true, cacheDir = cacdir);
    end
    engvars = dataGlider.parameterNames["eng"];
    scivars = dataGlider.parameterNames["sci"];

    # load engineering data from raw glider DBD files
    m_present_time = pyrow2jlcol(dataGlider.get("m_present_time"));
    m_lat = pyrow2jlcol(dataGlider.get("m_lat")); 
    m_lon = pyrow2jlcol(dataGlider.get("m_lon")); 
    m_pressure= pyrow2jlcol(dataGlider.get("m_pressure"));
    m_vacuum = pyrow2jlcol(dataGlider.get("m_vacuum"));
    m_battery = pyrow2jlcol(dataGlider.get("m_battery"));
    #m_leak_detect = pyrow2jlcol(dataGlider.get("m_leakdetect")); 
    m_veh_temp = pyrow2jlcol(dataGlider.get("m_veh_temp"));  
    m_roll = pyrow2jlcol(dataGlider.get("m_roll"));
    c_heading = pyrow2jlcol(dataGlider.get("c_heading")); 
    m_heading = pyrow2jlcol(dataGlider.get("m_heading")); 
    m_pitch = pyrow2jlcol(dataGlider.get("m_pitch"));
    m_battpos = pyrow2jlcol(dataGlider.get("m_battpos"));
    m_de_oil_vol = pyrow2jlcol(dataGlider.get("m_de_oil_vol"));
    m_ballast_pumped = pyrow2jlcol(dataGlider.get("m_ballast_pumped"));
    m_fin = pyrow2jlcol(dataGlider.get("m_fin"));
    m_depth_rate = pyrow2jlcol(dataGlider.get("m_depth_rate")); 
    m_altimeter_status = pyrow2jlcol(dataGlider.get("m_altimeter_status")); 
    m_raw_altitude = pyrow2jlcol(dataGlider.get("m_raw_altitude"));
    m_altimeter_voltage = pyrow2jlcol(dataGlider.get("m_altimeter_voltage")); 
    m_num_tot_inflections = pyrow2jlcol(dataGlider.get("m_tot_num_inflections")); 

    #m_ = pyrow2jlcol(dataGlider.get("m_")); 

    # load CTD data from raw glider EBD files
    #tctd, cond, temp, pres, m_de_oil_vol = dataGlider.get_CTD_sync("m_de_oil_vol");
    sci_m_present_time = pyrow2jlcol(dataGlider.get("sci_m_present_time"));
    sci_water_temp = pyrow2jlcol(dataGlider.get("sci_water_temp"));
    sci_water_cond = pyrow2jlcol(dataGlider.get("sci_water_cond"));
    sci_water_pressure = pyrow2jlcol(dataGlider.get("sci_water_pressure"));
    sci_flbbcd_chlor_units = pyrow2jlcol(dataGlider.get("sci_flbbcd_chlor_units"));
    sci_flbbcd_cdom_units = pyrow2jlcol(dataGlider.get("sci_flbbcd_cdom_units"));
    sci_flbbcd_bb_units = pyrow2jlcol(dataGlider.get("sci_flbbcd_bb_units"));
    sci_bsipar_par = pyrow2jlcol(dataGlider.get("sci_bsipar_par"));

    # calculate derived values from CTD data
    llat = Statistics.mean(m_lat[:,2]);
    llon = Statistics.mean(m_lon[:,2]);
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
    return ctdData
end

end