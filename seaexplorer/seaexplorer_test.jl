workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

using GLMakie, NCDatasets, NaNMath, Dates, Interpolations, MAT, GibbsSeaWater, Plots
import seaexplorer_functions: seaexplorer_MR_laur_load

if (@isdefined jm) == false
    display("Loading NORSE SEA064 data.")
    include("run_seaexplorer.jl")
end

pz = -20;
dz = 10;

pld1 = jmpld1d;
pld1t = datetime2unix.(pld1.t);

datadir = "/Users/gong/GitHub/jlglider/seaexplorer/data/"
mr1 = matopen(datadir * "jm_P.mat");
mr1_t = read(mr1, "unixt");
mr1_t[:,106] .= datetime2unix(DateTime(2022,10,30,1,0,0));
mr1_p = read(mr1, "pres");
mr1_eps1 = read(mr1, "eps1");
mr1_eps2 = read(mr1, "eps2");
mr1_z = gsw_z_from_p.(mr1_p, 70, 0, 0);

lon1f = linear_interpolation(pld1t, pld1.lon);
lat1f = linear_interpolation(pld1t, pld1.lat);

zlo, zhi = pz-dz, pz+dz;

lon1 = zeros(Float64, size(mr1_t,1), size(mr1_t,2));
lat1 = deepcopy(lon1);

eps1 = zeros(Float64, size(mr1_t,2));

for i = 1:size(mr1_p, 2)
    lon1[:,i] = lon1f(mr1_t[:,i]);
    lat1[:,i] = lat1f(mr1_t[:,i]);

    zind = findall((zlo .<= mr1_z[:,i] .<= zhi) .& (mr1_z[:,i] .< -1));
    eps1[i] = 10 .^ NaNMath.mean([NaNMath.mean(log10.(mr1_eps1[zind,i])); NaNMath.mean(log10.(mr1_eps2[zind,i]))]);
end

Plots.scatter(lon1[end,:], lat1[end,:]);