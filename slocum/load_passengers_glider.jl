using PyCall
using NaNMath, GibbsSeaWater, Dates, Interpolations
using GLMakie, ColorSchemes

dbdreader = pyimport("dbdreader");

# define function for converting an array from python row major to julia column major
function pyrow2jlcol(invar::Matrix{Float64})
    return reverse(rotr90(invar), dims = 2);
end

# specify valid data time period
datadir = "/Users/gong/Research/electa-20221103-passengers";
figoutdir = "/Users/gong/Research/electa-20221103-passengers/figures/";
t0 = DateTime("2022-10-01");
tN = DateTime("2022-12-01");
trange = datetime2unix.([t0; tN]);

# setup glider data loading using dbdreader
#dbdSylvia = dbdreader.DBD("/Users/gong/GitHub/sylvia-20180501-maracoos/data/DBD/01740010.DBD", cacheDir="/Users/gong/GitHub/sylvia-20180501-maracoos/data/cache/");
#dbdElecta = dbdreader.DBD(datadir * "/DBD/00970000.dbd", cacheDir = datadir * "/cache/");
#ebdElecta = dbdreader.DBD(datadir * "/EBD/00970000.ebd", cacheDir = datadir * "/cache/");

#dbdElecta = dbdreader.MultiDBD(datadir * "/DBD/0*.dbd", cacheDir = datadir * "/cache/");
#ebdElecta = dbdreader.MultiDBD(datadir * "/EBD/0*.ebd", cacheDir = datadir * "/cache/");

debdElecta = dbdreader.MultiDBD(pattern = datadir * "/[DE]BD/0*.[de]bd", complement_files = true, cacheDir = datadir * "/cache/");
dbdvars = debdElecta.parameterNames["eng"];
ebdvars = debdElecta.parameterNames["sci"];

# load CTD data from raw glider DBD & EBD files
#tctd, cond, temp, pres, m_de_oil_vol = debdElecta.get_CTD_sync("m_de_oil_vol");
sci_m_present_time = pyrow2jlcol(debdElecta.get("sci_m_present_time"));
sci_water_temp = pyrow2jlcol(debdElecta.get("sci_water_temp"));
sci_water_cond = pyrow2jlcol(debdElecta.get("sci_water_cond"));
sci_water_pressure = pyrow2jlcol(debdElecta.get("sci_water_pressure")); 

# QC loaded data by time and values, create interpolation function for T,C,P
tempind = findall((trange[1] .<= sci_water_temp[:,1] .<= trange[end]) .& (40.0 .>= sci_water_temp[:,2] .>= 0.0));
tempraw = sci_water_temp[tempind,:];
tempfunc = linear_interpolation(tempraw[:,1], tempraw[:,2], extrapolation_bc=Line());

condind = findall((trange[1] .<= sci_water_cond[:,1] .<= trange[end]) .& (100.0 .>= sci_water_cond[:,2] .>= 0.01));
condraw = sci_water_cond[condind,:];
condfunc = linear_interpolation(condraw[:,1], condraw[:,2], extrapolation_bc=Line());

presind = findall((trange[1] .<= sci_water_pressure[:,1] .<= trange[end]) .& (1000.0 .>= sci_water_pressure[:,2] .>= 0.01));
presraw = sci_water_pressure[presind,:];
presfunc = linear_interpolation(presraw[:,1], presraw[:,2], extrapolation_bc=Line());

# build a common timeline (not necessary if data is of good quality and CTD data all has the same length)
tall = unique(union(tempraw[:,1], condraw[:,1], presraw[:,1]));
tind = findall(trange[1] .<= tall .<= trange[end]);

t1 = floor(NaNMath.minimum(tall[tind]));
t2 = ceil(NaNMath.maximum(tall[tind]));
tt = collect(t1:1.0:t2);

# extract common time values for all CTD parameters
dtctd_fit = unix2datetime.(tall[tind]);
temp_fit = tempfunc(tall[tind]);
cond_fit = condfunc(tall[tind]);
pres_fit = presfunc(tall[tind]);
tctd_fit = deepcopy(tall[tind]);

gind = findall((30.0 .>= temp_fit .>= 15.0) .& (7.0 .>= cond_fit .>= 4.0) .& (50.0 .>= pres_fit .>= 0.0));
dtctd = dtctd_fit[gind];
tctd = tctd_fit[gind];
temp = temp_fit[gind];
cond = cond_fit[gind];
pres = pres_fit[gind];

si = sortperm(tctd);
tctd = tctd[si];
dtctd = dtctd[si];
temp = temp[si];
cond = cond[si];
pres = pres[si];

# calculate derived values from CTD data
gsw = GibbsSeaWater;
zz = gsw.gsw_z_from_p.(pres*10, 38.13, 0.0, 0.0);
salt = gsw.gsw_sp_from_c.(cond*10, temp, pres*10);
saltA = gsw.gsw_sa_from_sp.(salt, pres*10, -60.69, 38.13);
ctemp = gsw.gsw_ct_from_t.(saltA, temp, pres*10);
rho = gsw.gsw_rho.(saltA, ctemp, pres*10);
sigma0 = gsw.gsw_sigma0.(saltA, ctemp);
spice0 = gsw.gsw_spiciness0.(saltA, ctemp);
sndspd = gsw.gsw_sound_speed.(saltA, ctemp, pres*10);
