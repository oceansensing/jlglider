using Dates

mutable struct GliderGlobalAttributes
    title::String
    platform::String
    platform_vocabulary::String
    wmoid::String
    id::String
    naming_authority::String
    institution::String
    internal_mission_identifier::String
    geospatial_lat_min::Float64
    geospatial_lat_max::Float64
    geospatial_lon_min::Float64
    geospatial_lon_max::Float64
    geospatial_vertical_min::Float64
    geospatial_vertical_max::Float64
    time_coverage_start::String
    time_coverage_end::String
    site::String
    site_vocabulary::String
    program::Array{String}
    project::String
    network::String
    contributor_name::Array{String}
    contributor_email::Array{String}
    contributor_id::Array{String}
    contributor_role::Array{String}
    contributor_role_vocabulary::Array{String}
    agency::Array{String}
    agency_role::Array{String}
    agency_role_vocabulary::Array{String}
    agency_id::Array{String}
    agency_id_vocabulary::Array{String}
    uri::Array{String}
    doi::String
    rtqc_method::String
    rtqc_method_doi::String
    web_link::String
    comment::String
    data_created::String
    featureType::String
    Conventions::String
end

mutable struct GliderDims
    N_MEASUREMENTS::Int64
    N_PARAM::Int64
    N_SENSOR::Int64
end

mutable struct GliderLocation
    LATITUDE_GPS::Array{Float64}
    LONGITUDE_GPS::Array{Float64}
    TIME_GPS::Array{Union{Missing,DateTime}}

    LATITUDE::Array{Float64}
    LONGITUDE::Array{Float64}
    TIME::Array{Union{Missing,DateTime}}
end

mutable struct GliderInfo
    TRAJECTORY::String

    PLATFORM_TYPE::String
    PLATFORM_MODEL::String
    WMO_IDENIFIER::String

    PLATFORM_SERIAL_NUMBER::String
    PLATFORM_CODE::String
    PLATFORM_DEPTH_RATING::Int64

    ICES_CODE::String
    PLATFORM_MAKER::String

    DEPLOYMENT_DATE::String
    DEPLOYMENT_LATITUDE::String
    DEPLOYMENT_LONGITUDE::String

    FIELD_COMPARISON_REFERENCE::Array{String}

    GLIDER_FIRMWARE_VERSION::String
    LANDSTATION_VERSION::String
    BATTERY_TYPE::String
    BATTERY_PACK::String

    TELECOM_TYPE::String
    TRACKING_SYSTEM::String

    PHASE::String
    PHASE_QC::String
end

mutable struct GliderSensor
    SENSOR::String
    SENSOR_MAKER::String
    SENSOR_MODEL::String
    SENSOR_SERIAL_NUMBER::String
    SENSOR_CALIBRATION_DATE::String
end

mutable struct GliderParameter
    PARAMETER::String
    PARAMETER_SENSOR::String
    PARAMETER_UNITS::String
end

mutable struct GliderDataCTD
    TIME::Array{Any}
    TIME_CTD::Array{Any}
    LATITUDE::Array{Float64}
    LONGITUDE::Array{Float64}
    PRESSURE::Array{Float64}
    PRESSURE_QC::Array{Any}
    TEMPERATURE::Array{Float64}
    TEMPERATURE_QC::Array{Any}
    CONDUCTIVITY::Array{Float64}
    CONDUCTIVITY_QC::Array{Any}
    SALINITY::Array{Float64}
    SALINITY_QC::Array{Any}
end

mutable struct GliderDataDerived
    TIME::Array{Any}
    LATITUDE::Array{Float64}
    LONGITUDE::Array{Float64}
    PRESSURE::Array{Float64}
    PRESSURE_QC::Array{Any}
    Z::Array{Float64}
    ABSOLUTE_SALINITY::Array{Float64}
    CONSERVATIVE_TEMPERATURE::Array{Float64}
    DENSITY::Array{Float64}
    SIGMA0::Array{Float64}
    SPICE0::Array{Float64}
end

mutable struct GliderDataFLBBCD
    TIME::Array{Any}
    TIME_FLBBCD::Array{Any}
end
