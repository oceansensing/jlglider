include("MR_types.jl")
include("MR_io.jl")
include("moov.jl")

using Glob, MAT, JLD2, GibbsSeaWater
using .MR_types: MicroRiderRaw
#using .MR_io: MR_load_profile, MR_datasetup, MR_loadjld2


mrpath = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20231112-norse-janmayen-complete/mr_processing/jld2files/data_0003.jld2";
#mrp = MR_loadjld2(mrpath, 1);
mrr = load(mrpath);
mrp = mrr["mrprofile"];

#=
project = "NORSE"
mission = "JM"
year = 2023
profileid = [1; 2; 3; 4; 5];
glider = "SEA064"

project, mission, year
lat = 71.0

jld2dir, matdir, pdir = MR_datasetup(project::String, mission::String, year::Int)

moov("*.mat", pdir, matdir)

datafiles = Glob.glob("*.mat", matdir);
global mrdata = MicroRiderRaw[];

for i in profileid

    display(i)
    mrfile = matopen(datafiles[i]);

    fullPath = read(mrfile, "fullPath");
    fs_fast = read(mrfile, "fs_fast");
    fs_slow = read(mrfile, "fs_slow");
    header_version = read(mrfile, "header_version")
    t_slow = read(mrfile, "t_slow")
    t_fast = read(mrfile, "t_fast")
    setupfilestr = read(mrfile, "setupfilestr")
    cfgobj = read(mrfile, "cfgobj");
    header = read(mrfile, "header");
    filetime = read(mrfile, "filetime");
    date = read(mrfile, "date");
    time = read(mrfile, "time");
    Gnd = read(mrfile, "Gnd");
    Ax = read(mrfile, "Ax");
    Ay = read(mrfile, "Ay");
    T1 = read(mrfile, "T1");
    T1_dT1 = read(mrfile, "T1_dT1");
    T2 = read(mrfile, "T2");
    T2_dT2 = read(mrfile, "T2_dT2");
    sh1 = read(mrfile, "sh1");
    sh2 = read(mrfile, "sh2");
    P = read(mrfile, "P");
    P_dP = read(mrfile, "P_dP");
    PV = read(mrfile, "PV");
    V_Bat = read(mrfile, "V_Bat");
    Incl_Y = read(mrfile, "Incl_Y");
    Incl_X = read(mrfile, "Incl_X");
    Incl_T = read(mrfile, "Incl_T");
    odas_version = read(mrfile, "odas_version");
    vehicle_info = read(mrfile, "vehicle_info");
    t_fast_YD = read(mrfile, "t_fast_YD");
    t_slow_YD = read(mrfile, "t_slow_YD");
    Year = read(mrfile, "Year");
    Month = read(mrfile, "Month");
    Day = read(mrfile, "Day");
    Hour = read(mrfile, "Hour");
    Minute = read(mrfile, "Minute");
    Second = read(mrfile, "Second");
    Milli = read(mrfile, "Milli");
    T1_slow = read(mrfile, "T1_slow");
    T1_fast = read(mrfile, "T1_fast");
    T2_slow = read(mrfile, "T2_slow");
    T2_fast = read(mrfile, "T2_fast");
    P_slow = read(mrfile, "P_slow");
    P_fast = read(mrfile, "P_fast");
    z_slow = GibbsSeaWater.gsw_z_from_p.(P_slow, lat, 0.0, 0.0);
    z_fast = GibbsSeaWater.gsw_z_from_p.(P_fast, lat, 0.0, 0.0);
    temperature_fast = read(mrfile, "temperature_fast");
    W_slow = read(mrfile, "W_slow");
    W_fast = read(mrfile, "W_fast");
    speed_slow = read(mrfile, "speed_slow");
    speed_fast = read(mrfile, "speed_fast");
    gradT1 = read(mrfile, "gradT1");
    gradT2 = read(mrfile, "gradT2");
    input_parameters = read(mrfile, "input_parameters");
    params = read(mrfile, "params");

    mrprofile = MicroRiderRaw(fullPath, fs_fast, fs_slow, header_version, t_slow, t_fast, setupfilestr, cfgobj, header, filetime, date, time, Gnd, Ax, Ay, T1, T1_dT1, T2, T2_dT2, sh1, sh2, P, P_dP, PV, V_Bat, Incl_Y, Incl_X, Incl_T, odas_version, vehicle_info, t_fast_YD, t_slow_YD, Year, Month, Day, Hour, Minute, Second, Milli, T1_slow, T1_fast, T2_slow, T2_fast, P_slow, P_fast, z_slow, z_fast, temperature_fast, W_slow, W_fast, speed_slow, speed_fast, gradT1, gradT2, input_parameters, params);
    jldsave(jld2dir * basename.(datafiles[i])[1:end-4] * ".jld2", true; mrprofile);


    #mrr = load(mrpath)
end
=#