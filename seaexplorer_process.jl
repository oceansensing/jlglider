#include("seaexplorer_load_rt.jl")
using NaNMath
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt, cleanPress, cleanAD2CP, clean9999, cleanFLBBCDchl, cleanFLBBCDbb700, cleanFLBBCDcdom
import gsw_c2po: sigma0_from_t_sp

#t = cleanTime(sea064pld1d.t);
p = cleanPress(sea064pld1d.legato_pressure);
temp = cleanTemp(sea064pld1d.legato_temperature);
salt = cleanSalt(sea064pld1d.legato_salinity);
sigma0 = sigma0_from_t_sp(temp, salt, p, NaNMath.mean(sea064pld1d.lon), NaNMath.mean(sea064pld1d.lat))
mr_eps1 = cleanEPS(sea064pld1d.mr1000g_eps1);
mr_eps2 = cleanEPS(sea064pld1d.mr1000g_eps2);
mr_qc1 = clean9999(sea064pld1d.mr1000g_qc1);
mr_qc2 = clean9999(sea064pld1d.mr1000g_qc2);
mr_sh1_std = clean9999(sea064pld1d.mr1000g_sh1_std);
mr_sh2_std = clean9999(sea064pld1d.mr1000g_sh2_std);
mr_t1_avg = clean9999(sea064pld1d.mr1000g_t1_avg);
mr_t2_avg = clean9999(sea064pld1d.mr1000g_t2_avg);
ad2cp_Ueast = cleanAD2CP(sea064pld1d.ad2cp_Ueast);
ad2cp_Unorth = cleanAD2CP(sea064pld1d.ad2cp_Unorth);
ad2cp_Utot = cleanAD2CP(sea064pld1d.ad2cp_Utot);
ad2cp_Udir = cleanAD2CP(sea064pld1d.ad2cp_Udir);
ad2cp_qf = cleanAD2CP(sea064pld1d.ad2cp_qf);
chla = cleanFLBBCDchl(sea064pld1d.flbbcd_chl_scaled);
bb700 = cleanFLBBCDbb700(sea064pld1d.flbbcd_bb_700_scaled);
cdom = cleanFLBBCDcdom(sea064pld1d.flbbcd_cdom_scaled);
