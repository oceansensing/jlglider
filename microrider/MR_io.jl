module MR_io

using Glob, MAT, JLD2
import MR_types: MicroRiderRaw

function MR_datasetup(project, mission)
    if project == "NORSE"
        datadirJM = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/"
        datadirLBE = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/"
        if (mission == "JM") | (mission == "Jan Mayan")
            datadir = datadirJM;
        elseif (mission == "LBE") | (mission == "Lofoten Basin Eddy")
            datadir = datadirLBE;
        else
            display("Unknown mission. Exit 0.")
        end

        pdir = datadir * "mr1000g/";
        matdir = datadir * "mr1000g_processing/matfiles/";
        jld2dir = datadir * "mr1000g_processing/jld2files/";
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

function MR_mat2jld2(project::String, mission::String)
    jld2dir, matdir, pdir = MR_datasetup(project::String, mission::String)

    datafiles = Glob.glob("*.mat", matdir);

    global mrdata = MicroRiderRaw[];

    #mrzarr = Zarr.zgroup(datadir * "MR_data.zarr");

    Threads.@threads for i = 1:length(datafiles)
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
        temperature_fast = read(mrfile, "temperature_fast");
        W_slow = read(mrfile, "W_slow");
        W_fast = read(mrfile, "W_fast");
        speed_slow = read(mrfile, "speed_slow");
        speed_fast = read(mrfile, "speed_fast");
        gradT1 = read(mrfile, "gradT1");
        gradT2 = read(mrfile, "gradT2");
        input_parameters = read(mrfile, "input_parameters");
        params = read(mrfile, "params");

        mrprofile = MicroRiderRaw(fullPath, fs_fast, fs_slow, header_version, t_slow, t_fast, setupfilestr, cfgobj, header, filetime, date, time, Gnd, Ax, Ay, T1, T1_dT1, T2, T2_dT2, sh1, sh2, P, P_dP, PV, V_Bat, Incl_Y, Incl_X, Incl_T, odas_version, vehicle_info, t_fast_YD, t_slow_YD, Year, Month, Day, Hour, Minute, Second, Milli, T1_slow, T1_fast, T2_slow, T2_fast, P_slow, P_fast, temperature_fast, W_slow, W_fast, speed_slow, speed_fast, gradT1, gradT2, input_parameters, params);
        jldsave(jld2dir * datafiles[i][end-11:end-4] * ".jld2", true; mrprofile);
        #global mrdata = push!(mrdata, mrprofile);
        #close(mrfile)

        ## trying to use Zarr but not really working
        # Zarr.zopen(datadir * "MR_data.zarr", mode="w")
        # a = Zarr.zcreate(Float64, mrzarr, "fs_fast", 1)
        # a = fs_fast;
    end
end

function MR_loadjld2(jld2datafilepath)
    #datadirJM = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/"
    #datadirLBE = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/"
    #
    #datadir = datadirLBE;
    #pdir = datadir * "mr1000g/";
    #matdir = datadir * "mr1000g_processing/matfiles/";
    #jld2dir = datadir * "mr1000g_processing/jld2files/";

    mrr = MicroRiderRaw[];

    #if isdir(jld2dir)
    #    jldfiles = Glob.glob("*.jld2", jld2dir);
    #    if !isempty(jldfiles)
    #        i = 1
    #        display(i)
    #        data = jldopen(jldfiles[i], "r");
    #        mrr = data["mrprofile"];
    #    end #if
    #end #if

    data = jldopen(jld2datafilepath, "r");
    mrr = data["mrprofile"];
end

end #module