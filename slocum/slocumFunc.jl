# this module includes all the functions relevant for slocum glider data processing
# 2023-06-06  gong@vims.edu

module slocumFunc

using PyCall
using Glob, NaNMath, Statistics, GibbsSeaWater, Dates, Interpolations

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

function glider_presfunc(sci_water_pressure, trange)
    if isempty(sci_water_pressure) != true
        presind = findall((trange[1] .<= sci_water_pressure[1] .<= trange[end]) .& (1000.0 .>= sci_water_pressure[2] .>= 0.0));
        prestime = sci_water_pressure[1][presind]; 
        presraw = sci_water_pressure[2][presind];
        sortedpind = sortperm(prestime);
        presrawval = presraw[sortedpind];
        prestime = prestime[sortedpind];
        presfunc = linear_interpolation(prestime, presrawval, extrapolation_bc=Line());
        #presdtime = unix2datetime.(prestime);
        #prespres = presraw[:,2];
        #presz = gsw.gsw_z_from_p.(prespres*10, llat, 0.0, 0.0); 
    else
        presfunc = []; # need a function that returns [] when called with any input
        prestime = [];
        presraw = [];
    end
    return presfunc, prestime, presraw
end

# this function takes glider variable from dbdreader as input and format it for analysis (trimming, outlier removal, calculate)

function glider_var_load(varin, varrange)
    if isempty(varin) != true
        varind = findall(varrange[1] .<= varin .<= varrange[end]); 
        varout = glidervar[varind];
    else
        varind = [];
        varout = [];
    end
    return varind, varout;
end

function glider_ctd_load(t::Array{Float64}, pres::Array{Float64}, temp::Array{Float64}, cond::Array{Float64}, trange::Array{Float64})
    gind = findall((trange[1] .<= t .<= trange[end]) .& (0.0 .<= pres .<= 1000.0) .& (0.1 .<= temp .<= 40.0) .& (0.01 .<= cond .<= 100.0)); 
    return t[gind], pres[gind], temp[gind], cond[gind];
end

function glider_var_load(t::Array{Float64}, pres::Array{Float64}, glidervar::Array{Float64}, trange::Array{Float64}, varlim::Array{Float64})
    if isempty(glidervar) != true
        varind = findall((trange[1] .<= t .<= trange[end]) .& (0.0 .<= pres .<= 1000.0));
        t = t[varind];
        varout = glidervar[varind];
        presout = pres[varind];
        varbadind = findall((varlim[1] .>= varout) .| (varout .>= varlim[end]));
        varout[varbadind] .= NaN;
    else
        varind = [];
        t = [];
        varout = [];
        presout = [];
    end
    return t, presout, varout, varind;
end

function glider_var_load(glidervar, trange, varlim, p, lat)
    gsw = GibbsSeaWater;
    if isempty(glidervar) != true
        presfunc, prestime, presraw = glider_presfunc(p, trange);
        varind = findall((trange[1] .<= glidervar[1] .<= trange[end]) .& (varlim[1] .<= glidervar[2] .<= varlim[end])); 
        vartime = glidervar[1][varind];
        varraw = glidervar[2][varind];
        sortedvarind = sortperm(vartime);
        varrawval = varraw[sortedvarind];
        vartime = vartime[sortedvarind];
        varfunc = linear_interpolation(vartime, varrawval, extrapolation_bc=Line()); 
        vardtime = unix2datetime.(vartime);
        varpres = presfunc(vartime);
        varz = gsw.gsw_z_from_p.(varpres*10, lat, 0.0, 0.0);  
    else
        varraw = [];
        vartime = [];
        varpres = [];
        varz = [];
    end
    return vartime, varpres, varraw, varz
end

function datetick(unix_t)
    x = unix_t
    xdt = unix2datetime.(unix_t); 
    yday = Dates.dayofyear.(xdt);
    uyday = unique(yday);
    hour = Dates.Hour.(xdt);
    minute = Dates.Minute.(xdt);
    #df = DateFormat("y-m-d");

    #tickind = Vector{Int64}(undef, length(uyday));
    tickind = [];
    for i = 1:length(uyday)
        t0 = findall((yday .== uyday[i]) .& (hour .== Hour(0)));
        if isempty(t0) == false
            push!(tickind, t0[1]);
        end
    end
    xtick = x[tickind];
    #xticklabel = string.(Dates.Date.(xdt[tickind]));
    xticklabel = [x[6:10] for x in string.(xdt[tickind])];
    return xdt, xtick, xticklabel    
end

end #slocumFunc module