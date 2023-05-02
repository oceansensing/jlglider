# load roughly processed MR1000G data by Laur Ferris for plotting
# gong@vims.edu 2023-04-28

using MAT, GibbsSeaWater, Interpolations, NaNMath

#function freya_MR_laur_load(pld1, pld2, pz, dz)


datadir = "/Users/gong/GitHub/jlglider/seaexplorer/data/"
glider = matread(datadir * "freya.mat");

gtime = glider["freya"]["unixt"][:];
gpres = glider["freya"]["sci_water_pressure"][:];
gtemp = glider["freya"]["sci_water_temp"][:];
gcond = glider["freya"]["sci_water_cond"][:];
gpres_tp = glider["freya"]["sci_lsltp_pressure"][:];
gtime_tp = glider["freya"]["sci_lsltp_timestamp"][:];
geps1_tp = glider["freya"]["sci_lsltp_e1"][:];
geps2_tp = glider["freya"]["sci_lsltp_e2"][:];

gnavlon = glider["freyanav"]["lon"][:];
gnavlat = glider["freyanav"]["lat"][:];
gnavtime = glider["freyanav"]["unixt"][:];

glonf = linear_interpolation(gnavtime, gnavlon, extrapolation_bc=Line());
glatf = linear_interpolation(gnavtime, gnavlat, extrapolation_bc=Line());

#zlo, zhi = pz-dz, pz+dz;

glon = glonf(gtime[:]);
glat = glatf(gtime[:]);

gz = gsw_z_from_p.(gpres, 70, 0, 0);
gsalt = gsw_sp_from_c.(gcond, gtemp, gpres);
gsaltA = gsw_sa_from_sp.(gsalt, gpres, glon, glat);
gctemp = gsw_ct_from_t.(gsaltA, gtemp, gpres);
gsigma0 = gsw_sigma0.(gsaltA, gctemp);
gspice0 = gsw_spiciness0.(gsaltA, gctemp);
gsndspd = gsw_sound_speed.(gsaltA, gctemp, gpres);



#end
    