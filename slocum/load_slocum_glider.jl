# this script loads Slocum glider data as configured for MARACOOS W missions using dbdreader
# gong@vims.edu 2023-03-26: adopted from the PASSENGERS version - added sorting of raw data by time and plotting of chla data
#

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations

dbdreader = pyimport("dbdreader");
gsw = GibbsSeaWater;

mutable struct engStruct
    t::Array{DateTime}
    p::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
end
mutable struct sciStruct
    t::Array{DateTime}
    p::Array{Float64}
    z::Array{Float64}
    lon::Array{Float64}
    lat::Array{Float64}
    temp::Array{Float64}
    cond::Array{Float64}
    salt::Array{Float64}
    ctemp::Array{Float64}
    saltA::Array{Float64}
    sigma0::Array{Float64}
    spice0::Array{Float64}
    sndspd::Array{Float64}
    castnum::Int
    qflag::Int
end

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
#dbdSylvia = dbdreader.DBD("/Users/gong/GitHub/sylvia-20180501-maracoos/data/DBD/01740010.DBD", cacheDir="/Users/gong/GitHub/sylvia-20180501-maracoos/data/cache/");
#dbdElecta = dbdreader.DBD(datadir * "/DBD/00970000.dbd", cacheDir = datadir * "/cache/");
#ebdElecta = dbdreader.DBD(datadir * "/EBD/00970000.ebd", cacheDir = datadir * "/cache/");

#dbdElecta = dbdreader.MultiDBD(datadir * "/DBD/0*.dbd", cacheDir = datadir * "/cache/");
#ebdElecta = dbdreader.MultiDBD(datadir * "/EBD/0*.ebd", cacheDir = datadir * "/cache/");

#debdElecta = dbdreader.MultiDBD(pattern = datadir * "/[DE]BD/0*.[de]bd", complement_files = true, cacheDir = datadir * "/cache/");
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

# QC loaded data by time and values, create interpolation function for T,C,P
tempind = findall((trange[1] .<= sci_water_temp[:,1] .<= trange[end]) .& (40.0 .>= sci_water_temp[:,2] .> 0.0));
tempraw = sci_water_temp[tempind,:];
sortedtind = sortperm(tempraw[:,1]);
tempraw = tempraw[sortedtind,:];
tempfunc = linear_interpolation(tempraw[:,1], tempraw[:,2], extrapolation_bc=Line());

condind = findall((trange[1] .<= sci_water_cond[:,1] .<= trange[end]) .& (100.0 .>= sci_water_cond[:,2] .>= 0.01));
condraw = sci_water_cond[condind,:];
sortedcind = sortperm(condraw[:,1]);
condraw = condraw[sortedcind,:];
condfunc = linear_interpolation(condraw[:,1], condraw[:,2], extrapolation_bc=Line());

presind = findall((trange[1] .<= sci_water_pressure[:,1] .<= trange[end]) .& (1000.0 .>= sci_water_pressure[:,2] .>= 0.01));
presraw = sci_water_pressure[presind,:];
sortedpind = sortperm(presraw[:,1]);
presraw = presraw[sortedpind,:];
presfunc = linear_interpolation(presraw[:,1], presraw[:,2], extrapolation_bc=Line());

chlaind = findall((trange[1] .<= sci_flbbcd_chlor_units[:,1] .<= trange[end]) .& (2.0 .>= sci_flbbcd_chlor_units[:,2] .>= -0.1)); 
chlaraw = sci_flbbcd_chlor_units[chlaind,:];
sortedchlaind = sortperm(chlaraw[:,1]);
chlaraw = chlaraw[sortedchlaind,:];
chlafunc = linear_interpolation(chlaraw[:,1], chlaraw[:,2], extrapolation_bc=Line()); 
chlatime = chlaraw[:,1];
chladtime = unix2datetime.(chlatime);
chlapres = presfunc(chlatime);
chlaz = gsw.gsw_z_from_p.(chlapres*10, llat, 0.0, 0.0); 

cdomind = findall((trange[1] .<= sci_flbbcd_cdom_units[:,1] .<= trange[end]) .& (5.0 .>= sci_flbbcd_cdom_units[:,2] .>= -5.0)); 
cdomraw = sci_flbbcd_cdom_units[cdomind,:];
sortedcdomind = sortperm(cdomraw[:,1]);
cdomraw = cdomraw[sortedcdomind,:];
cdomfunc = linear_interpolation(cdomraw[:,1], cdomraw[:,2], extrapolation_bc=Line()); 
cdomtime = cdomraw[:,1];
cdomdtime = unix2datetime.(cdomtime);
cdompres = presfunc(cdomtime);
cdomz = gsw.gsw_z_from_p.(cdompres*10, llat, 0.0, 0.0); 

bb700ind = findall((trange[1] .<= sci_flbbcd_bb_units[:,1] .<= trange[end]) .& (0.008 .>= sci_flbbcd_bb_units[:,2] .>= 0.0)); 
bb700raw = sci_flbbcd_bb_units[bb700ind,:];
sortedbb700ind = sortperm(bb700raw[:,1]);
bb700raw = bb700raw[sortedbb700ind,:];
bb700func = linear_interpolation(bb700raw[:,1], bb700raw[:,2], extrapolation_bc=Line()); 
bb700time = bb700raw[:,1];
bb700dtime = unix2datetime.(bb700time);
bb700pres = presfunc(bb700time);
bb700z = gsw.gsw_z_from_p.(bb700pres*10, llat, 0.0, 0.0); 

bparind = findall((trange[1] .<= sci_bsipar_par[:,1] .<= trange[end]) .& (6000.0 .>= sci_bsipar_par[:,2] .>= 0)); 
bparraw = sci_bsipar_par[bparind,:];
sortedbparind = sortperm(bparraw[:,1]);
bparraw = bparraw[sortedbparind,:];
bparfunc = linear_interpolation(bparraw[:,1], bparraw[:,2], extrapolation_bc=Line()); 
bpartime = bparraw[:,1];
bpardtime = unix2datetime.(bpartime);
bparpres = presfunc(bpartime);
bparz = gsw.gsw_z_from_p.(bparpres*10, llat, 0.0, 0.0); 

# find common glider values
tctd = unique(intersect(presraw[:,1], tempraw[:,1], condraw[:,1]));
tctdT = intersectalajulia2(tctd, tempraw[:,1])[3];
tctdC = intersectalajulia2(tctd, condraw[:,1])[3];
tctdP = intersectalajulia2(tctd, presraw[:,1])[3];

tpuck = chlaraw[:,1];

# build a common timeline (not necessary if data is of good quality and CTD data all has the same length)
tall = unique(union(tempraw[:,1], condraw[:,1], presraw[:,1]));
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

gind = findall((30.0 .>= temp_fit .>= 0.0) .& (7.0 .>= cond_fit .>= 3.0) .& (50.0 .>= pres_fit .>= 0.0) .& (20.0 .>= chla_fit .>= 0));
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

# values fitted to a 1 second time grid
zzf = gsw.gsw_z_from_p.(presf*10, llat, 0.0, 0.0);
saltf = gsw.gsw_sp_from_c.(condf*10, tempf, presf*10);
saltAf = gsw.gsw_sa_from_sp.(saltf, presf*10, llon, llat);
ctempf = gsw.gsw_ct_from_t.(saltAf, tempf, presf*10);
rhof = gsw.gsw_rho.(saltAf, ctempf, presf*10);
sigma0f = gsw.gsw_sigma0.(saltAf, ctempf);
spice0f = gsw.gsw_spiciness0.(saltAf, ctempf);
sndspdf = gsw.gsw_sound_speed.(saltAf, ctempf, presf*10);

engData = engStruct[];
sciData = sciStruct[];