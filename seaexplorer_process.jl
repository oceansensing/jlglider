#include("seaexplorer_load_rt.jl")
import seaexplorer_functions: cleanEPS, cleanTemp, cleanSalt

using NaNMath, Dates

eps1 = cleanEPS(sea064pld1d.mr1000g_eps1);
eps2 = cleanEPS(sea064pld1d.mr1000g_eps2);
temp = cleanTemp(sea064pld1d.legato_temperature);
salt = cleanSalt(sea064pld1d.legato_salinity);
