#include("seaexplorer_load_rt.jl")
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt, cleanFLBBCDchl, cleanFLBBCDbb700, cleanFLBBCDcdom

eps1 = cleanEPS(sea064pld1d.mr1000g_eps1);
eps2 = cleanEPS(sea064pld1d.mr1000g_eps2);
temp = cleanTemp(sea064pld1d.legato_temperature);
salt = cleanSalt(sea064pld1d.legato_salinity);
chla = cleanFLBBCDchl(sea064pld1d.flbbcd_chl_scaled);
bb700 = cleanFLBBCDbb700(sea064pld1d.flbbcd_bb_700_scaled);
cdom = cleanFLBBCDcdom(sea064pld1d.flbbcd_cdom_scaled);
