module MR_io

include("MR_types.jl")
include("moov.jl")

using Glob, MAT, JLD2, GibbsSeaWater
using .MR_types: MicroRiderRaw

function MR_datasetup(project, mission, year)
    if project == "NORSE"
        datadirJM2022 = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/"
        datadirLBE2022 = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/"
        datadirJM2023 = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20231112-norse-janmayen-complete/"
        
        if year == 2022
            if (mission == "JM") | (mission == "Jan Mayan")
                datadir = datadirJM2022;
            elseif (mission == "LBE") | (mission == "Lofoten Basin Eddy")
                datadir = datadirLBE2022;
            else
                display("Unknown mission. Exit 0.")
            end

            pdir = datadir * "mr1000g/";
            matdir = datadir * "mr1000g_processing/matfiles/";
            jld2dir = datadir * "mr1000g_processing/jld2files/";
        elseif year == 2023
            datadir = datadirJM2023;
            pdir = datadir * "mr/";
            matdir = datadir * "mr_processing/matfiles/";
            jld2dir = datadir * "mr_processing/jld2files/";
        else
            display("Unknown year. Exit 0.")
        end
    end

    if !isdir(matdir)
        mkpath(matdir)
        println("matdir directory created.")
    end

    if !isdir(jld2dir)
        mkpath(jld2dir)
        println("jld2dir directory created.")
    end
    return jld2dir, matdir, pdir
end

function MR_datasetup(project, mission)
    MR_datasetup(project, mission, 2022)
end

# DO NOT use this function! use the MR_mat2jld2.jl script instead. there is a type import/export issue with jldsave. 2024-04-04 DG
function MR_mat2jld2(project::String, mission::String, year::Int, lat::Float64)
    jld2dir, matdir, pdir = MR_datasetup(project::String, mission::String, year::Int)

    moov("*.mat", pdir, matdir)

    datafiles = Glob.glob("*.mat", matdir);
    global mrdata = MicroRiderRaw[];

    Threads.@threads for i = 1:length(datafiles)
    #for i = 11:15

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
    end
end

function MR_loadjld2(jld2datafilepath::String; loadflag = 1)
    #if (@isdefined loadflag) == 0
    #    loadflag = 1;
    #end

    if loadflag == 0 # for some reason, running jlopen more than once lead to type issues when using push!(), 20240404 DG
        mrr = jldopen(jld2datafilepath, "r");
        return mrr["mrprofile"];
    else # this is the preferred way to load the MR .jld2 files, 20240404 DG
        mrr = load(jld2datafilepath);
        return mrr["mrprofile"];
    end
end

# load specific profiles based on what user desires
function MR_load_profile(project::String, mission::String, year::Int, profilename::String)
    display(profilename)
    jld2dir, matdir, pdir = MR_datasetup(project, mission, year);
    datafiles = Glob.glob("*.jld2", jld2dir);
    filenames = String[];
    for i = 1:length(datafiles)
        push!(filenames, basename.(datafiles[i])[1:end-5]);
    end
    iprofile = findall(filenames .== profilename);
    if !isempty(iprofile)
        MRprofile = MR_loadjld2(datafiles[iprofile[1]]; loadflag = 1);
    else
        display("Profile " * profilename * " not found, there are " * string(length(filenames)) * " available, exit 0.")
        MRprofile = 0;
    end
    return MRprofile;
end

function MR_load_profile(project::String, mission::String, year::Int, profileid::Int)
    if year == 2022
        datafilestr = "dat_" * string(profileid, base = 10, pad = 4);
    elseif year == 2023
        datafilestr = "data_" * string(profileid, base = 10, pad = 4);
    else
        display("Unknown year. Exit 0.")
    end
    MRprofile = MR_load_profile(project, mission, year, datafilestr);
    return MRprofile;
end

end #module