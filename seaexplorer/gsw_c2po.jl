module gsw_c2po

using GibbsSeaWater

function rho_from_t_sp(t, sp, p, lon, lat)
    saltA = gsw_sa_from_sp.(sp, p, lon, lat);
    ctemp = gsw_ct_from_t.(saltA, t, p);
    sigma0 = gsw_sigma0.(saltA, ctemp);
    rho = gsw_rho.(saltA, ctemp, p);
    return rho
end

function sigma0_from_t_sp(t, sp, p, lon, lat)
    saltA = gsw_sa_from_sp.(sp, p, lon, lat);
    ctemp = gsw_ct_from_t.(saltA, t, p);
    sigma0 = gsw_sigma0.(saltA, ctemp);
    return sigma0;
end

function spice0_from_t_sp(t, sp, p, lon, lat)
    saltA = gsw_sa_from_sp.(sp, p, lon, lat);
    ctemp = gsw_ct_from_t.(saltA, t, p);
    spice0 = gsw_spiciness0.(saltA, ctemp);
    return spice0;
end

function N2_from_t_sp(t, sp, p, lon, lat)
    saltA = gsw_sa_from_sp.(sp, p, lon, lat);
    ctemp = gsw_ct_from_t.(saltA, t, p);
    n2, pmid = gsw_nsquared(saltA, ctemp, p, lat);
end

end