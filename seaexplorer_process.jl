#include("seaexplorer_load_rt.jl")
using NaNMath
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt, cleanPress, cleanFLBBCDchl, cleanFLBBCDbb700, cleanFLBBCDcdom
import gsw_c2po: sigma0_from_t_sp

eps1 = cleanEPS(sea064pld1d.mr1000g_eps1);
eps2 = cleanEPS(sea064pld1d.mr1000g_eps2);
temp = cleanTemp(sea064pld1d.legato_temperature);
salt = cleanSalt(sea064pld1d.legato_salinity);
p = cleanPress(sea064pld1d.legato_pressure);
sigma0 = sigma0_from_t_sp(temp, salt, p, NaNMath.mean(sea064pld1d.lon), NaNMath.mean(sea064pld1d.lat))
chla = cleanFLBBCDchl(sea064pld1d.flbbcd_chl_scaled);
bb700 = cleanFLBBCDbb700(sea064pld1d.flbbcd_bb_700_scaled);
cdom = cleanFLBBCDcdom(sea064pld1d.flbbcd_cdom_scaled);
