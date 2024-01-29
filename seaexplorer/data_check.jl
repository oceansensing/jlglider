using Glob, DataFrames, CSV, Dates, Missings, NaNMath, Interpolations, Statistics
using NCDatasets
using GibbsSeaWater, MAT
using Plots

jm23_gli = CSV.read("/Users/rbbourdon/Desktop/NORSE/ad2cp_processed/M48_glider_data_processed.csv", DataFrame)
jm23_sci = CSV.read("/Users/rbbourdon/Desktop/NORSE/ad2cp_processed/M48_science_data_processed.csv", DataFrame)

jm23_sci.density = gsw.density.rho(jm23_sci.SA, jm23_sci.CT, jm23_sci.depth)

function microrider_check(data::DataFrame, turbulence::String, dive_num::Int)
    gdf = groupby(data, :diveNum)
    dives = unique(data.diveNum)

    yo = gdf[dive_num]
    avg_turb = mean(yo[turbulence])
    std_turb = std(yo[turbulence])

    println("The mean for $turbulence is $avg_turb w/ a std of +/- $std_turb")

    plot = Plots.scatter(yo[turbulence], yo.depth, label = turbulence, xlabel = turbulence, ylabel = "Depth", yflip = true, title = "Turbulence vs Depth for Dive $dive_num")

    # Show the plot
    display(plot)
end 

function ctd_data_check(data::DataFrame, dive_num::Int)
    gdf = groupby(data, :diveNum)
    dives = unique(data.diveNum)

    yo = gdf[dive_num]
    sigma0 = gsw.density.sigma0(yo.SA, yo.CT)
    y = yo.depth
    # Create a subplot layout with 1 row and 3 columns
    layout = @layout [a b c]

    # Create the first subplot
    p1 = scatter(yo.CT, y, label = "Conserative Temperature", c = "red", xlabel = "Conserative Temperature", ylabel = "Depth", yflip = true)

    # Create the second subplot
    p2 = plot(yo.SA, y, label= "Absolute Salinity", c = "blue", xlabel = "Absolute Salinity", yflip = true)

    # Create the third subplot
    p3 = plot(sigma0, y, label = "Sigma", c = "green", xlabel = "Desnity (Sigma)", yflip = true)

    # Combine the subplots into a single plot
    plot = plot(p1, p2, p3, layout = layout)
    display(plot)
end

function tu(data::DataFrame, window::Int)
    global turner
    turner_list = []
    gdf = groupby(data, :diveNum)
    dives = unique(data.diveNum)

    for dive in dives
        yo = gdf[dive]
        dz = median(diff(yo.depth))
        window_n = Int(window // dz)
        if window_n < 0
            window_n = -1 * window_n
        end
        CT = movmean(yo.CT, window_n)
        SA = movmean(yo.SA, window_n)
        Tu, R_rho, p_mid = gsw.stability.Turner_Rsubrho(SA, CT, yo.depth)
        yo1 = yo[2:end, :]
        nanind = findall(!isnan.(Tu))
        Tu = Tu[nanind]
        R_rho = R_rho[nanind]
        p_mid = p_mid[nanind]
        time = yo1.date_float[nanind]
        d = DataFrame(Tu = Tu, R_rho = R_rho, p_mid = p_mid, Time = time)
        turn = d
        dt_tot = push!(turner_list, turn)
        turner = vcat(turner_list...)
    end
end

function n2(data::DataFrame, window::Int)
    global n2
    n2_list = []
    gdf = groupby(data, :diveNum)
    dives = unique(data.diveNum)

    for dive in dives
        yo = gdf[dive]
        dz = median(diff(yo.depth))
        window_n = Int(window // dz)
        if window_n < 0
            window_n = -1 * window_n
        end
        CT = movmean(yo.CT, window_n)
        SA = movmean(yo.SA, window_n)
        N2, p_mid = gsw.stability.Nsquared(SA, CT, yo.depth, lat = yo.Lat)
        yo1 = yo[2:end, :]
        nanind = findall(!isnan.(N2))
        N2 = N2[nanind]
        p_mid = p_mid[nanind]
        time = yo1.date_float[nanind]
        d = DataFrame(N2 = N2, p_mid = p_mid, Time = time)
        n2 = d
        dt_tot = push!(n2_list, n2)
        n2 = vcat(n2_list...)
    end
end
