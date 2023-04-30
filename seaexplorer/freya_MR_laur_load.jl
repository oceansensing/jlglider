# load roughly processed MR1000G data by Laur Ferris for plotting
# gong@vims.edu 2023-04-28

using MAT, GibbsSeaWater, Interpolations, NaNMath

function freya_MR_laur_load(pld1, pld2, pz, dz)

datadir = "/Users/gong/GitHub/jlglider/seaexplorer/data/"
freya = matopen(datadir * "freya.mat");

mr1_t = read(mr1, "unixt");
mr1_p = read(mr1, "pres");
mr1_eps1 = read(mr1, "eps1");
mr1_eps2 = read(mr1, "eps2");

mr2_t = read(mr2, "unixt");
mr2_p = read(mr2, "pres");
mr2_eps1 = read(mr2, "eps1");
mr2_eps2 = read(mr2, "eps2");

mr1_z = gsw_z_from_p.(mr1_p, 70, 0, 0);
mr2_z = gsw_z_from_p.(mr2_p, 70, 0, 0);

lon1f = linear_interpolation(pld1t, pld1.lon);
lat1f = linear_interpolation(pld1t, pld1.lat);
lon2f = linear_interpolation(pld2t, pld2.lon);
lat2f = linear_interpolation(pld2t, pld2.lat); #  extrapolation_bc=Line()

zlo, zhi = pz-dz, pz+dz;

lon1 = zeros(Float64, size(mr1_t,1), size(mr1_t,2));
lat1 = deepcopy(lon1);
lon2 = zeros(Float64, size(mr2_t,1), size(mr2_t,2));
lat2 = deepcopy(lon2);

eps1 = zeros(Float64, size(mr1_t,2));
eps2 = zeros(Float64, size(mr2_t,2));

for i = 1:size(mr1_p, 2)
    lon1[:,i] = lon1f(mr1_t[:,i]);
    lat1[:,i] = lat1f(mr1_t[:,i]);

    zind = findall(zlo .<= mr1_z[:,i] .<= zhi);
    eps1[i] = 10 .^ NaNMath.mean(log10.(mr1_eps1[zind,i]));
end

for i = 1:size(mr2_p, 2)
    lon2[:,i] = lon2f(mr2_t[:,i]);
    lat2[:,i] = lat2f(mr2_t[:,i]);

    zind = findall(zlo .<= mr2_z[:,i] .<= zhi);
    eps2[i] = 10 .^ NaNMath.mean(log10.(mr2_eps1[zind,i]));
end

return lon1, lat1, lon2, lat2, eps1, eps2
end
    