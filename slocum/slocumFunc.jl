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

function glider_presfunc(sci_water_pressure, trange::Array{Float64})
    if isempty(sci_water_pressure) != true
        presind = findall((trange[1] .<= sci_water_pressure[1] .<= trange[end]) .& (300.0 .>= sci_water_pressure[2] .>= 0.0));
        prestime = sci_water_pressure[1][presind]; 
        presraw = sci_water_pressure[2][presind];
        sortedpind = sortperm(prestime);
        presrawval = presraw[sortedpind];
        prestime = prestime[sortedpind];
        presfunc = linear_interpolation(prestime, presrawval, extrapolation_bc=Line());
        #BSpline(Linear())
        #presfunc = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate(presrawval, BSpline(Linear())), prestime));
        #presdtime = unix2datetime.(prestime);
        #prespres = presraw[:,2];
        #presz = gsw.gsw_z_from_p.(prespres*10, llat, 0.0, 0.0); 
    else
        presfunc = []; # need a function that returns [] when called with any input
        prestime = [];
        presrawval = [];
    end
    return presfunc, prestime, presrawval
end

function glider_presfunc(sci_water_pressure, tbound::Float64)
    if isempty(sci_water_pressure) != true
        presind = findall((median(sci_water_pressure[1]) - tbound .<= sci_water_pressure[1] .<= median(sci_water_pressure[1]) + tbound) .& (300.0 .>= sci_water_pressure[2] .>= 0.0));
        prestime = sci_water_pressure[1][presind]; 
        presraw = sci_water_pressure[2][presind];
        sortedpind = sortperm(prestime);
        presrawval = presraw[sortedpind];
        prestime = prestime[sortedpind];
        presfunc = linear_interpolation(prestime, presrawval, extrapolation_bc=Line());
        #BSpline(Linear())
        #presfunc = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate(presrawval, BSpline(Linear())), prestime));
        #presdtime = unix2datetime.(prestime);
        #prespres = presraw[:,2];
        #presz = gsw.gsw_z_from_p.(prespres*10, llat, 0.0, 0.0); 
    else
        presfunc = []; # need a function that returns [] when called with any input
        prestime = [];
        presrawval = [];
    end
    return presfunc, prestime, presrawval
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

#=
function glider_ctd_qc(t::Array{Float64}, pres::Array{Float64}, temp::Array{Float64}, cond::Array{Float64}, trange::Array{Float64})
    gind = findall((trange[1] .<= t .<= trange[end]) .& (0.0 .<= pres .<= 1000.0) .& (17.5 .<= temp .<= 40.0) .& (4.0 .<= cond .<= 6.0)); 
    return t[gind], pres[gind], temp[gind], cond[gind];
end
=#

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

function glider_var_load(glidervar, trange::Array{Float64}, varlim, sci_water_pressure, lat)
    gsw = GibbsSeaWater;
    if isempty(glidervar) != true
        presfunc, prestime, presraw = glider_presfunc(sci_water_pressure, trange);
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
        varfunc = [];
    end
    return varfunc, vartime, varpres, varrawval, varz
end

function glider_var_load(glidervar, tbound::Float64, varlim, presfunc, lat)
    gsw = GibbsSeaWater;
    if isempty(glidervar) != true
        #presfunc, prestime, presraw = glider_presfunc(p, trange);
        varind = findall((median(glidervar[1]) - tbound .<= glidervar[1] .<= median(glidervar[1]) + tbound) .& (varlim[1] .<= glidervar[2] .<= varlim[end])); 
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
        varrawval = [];
        vartime = [];
        varpres = [];
        varz = [];
        varfunc = [];
    end
    return varfunc, vartime, varrawval, varpres, varz
end


# DG 2024-09-06, adopted from seaexplorer
function datetime2yearday(xdt::DateTime)
    year = Dates.year(xdt);
    yday = Dates.dayofyear(xdt);
    seconds_in_day = 86400;
    ydayfrac = yday + (Dates.hour(xdt) * 3600 .+ Dates.minute(xdt) * 60 .+ Dates.second(xdt)) ./ seconds_in_day;
    return ydayfrac;
end

# DG 2024-09-06, ChatGPT 4o
function yearday2datetime(yyyy::Int, yearday::Float64)
    # Separate integer part (day) and fractional part (time of day)
    int_day = Int(floor(yearday))  # Get the integer part of the day
    fractional_day = yearday - int_day  # Get the fractional part of the day

    # Convert fractional day into hours, minutes, and seconds
    total_seconds = Int(round(fractional_day * 86400))  # 86400 seconds in a day
    hours = div(total_seconds, 3600)  # Number of full hours
    minutes = div(total_seconds % 3600, 60)  # Number of full minutes
    seconds = total_seconds % 60  # Remaining seconds

    # Start from January 1st of the given year and add the number of days
    dateint = Date(yyyy, 1, 1) + Day(int_day - 1)  # Subtract 1 because days start from 1

    # Create the final DateTime with the computed hours, minutes, and seconds
    return DateTime(year(dateint), month(dateint), day(dateint), hours, minutes, seconds)
    #return DateTime(date) + Time(hours, minutes, seconds)
end

end #slocumFunc module