# Julia implementation of Laur Ferris' glider_debd_processing.m
# Donglai Gong 2020-10-05

using Glob, Dates

mutable struct SensorList
    transmitted::Array{Bool}
    sensnum::Array{Int}
    indexnum::Array{Int}
    nbytes::Array{Int}
    name::Array{String}
    unit::Array{String}
end

# specifying data directory
dbd_dir = "/path/to/glider/data/DBD/";
ebd_dir = "/path/to/glider/data/EBD/";

sbd_dir = "/path/to/glider/data/SBD/";
tbd_dir = "/path/to/glider/data/TBD/";

# specifying cache directory
cac_dir = "/Users/gong/GitHub/data_processing/cache/";
#cache_path_dbd = "/Users/gong/GitHub/data_processing/cache/be93efad.cac";
#cache_path_ebd = "/Users/gong/GitHub/data_processing/cache/2804f449.cac";

# obtaining the paths to all cache and glider data files
cac_path1 = Glob.glob("*.cac", cac_dir);
cac_path2 = Glob.glob("*.CAC", cac_dir);
cac_path = [cac_path1; cac_path2];

dbd_path1 = Glob.glob("*.dbd", dbd_dir);
dbd_path2 = Glob.glob("*.DBD", dbd_dir);
dbd_path = [dbd_path1; dbd_path2];

ebd_path1 = Glob.glob("*.ebd", ebd_dir);
ebd_path2 = Glob.glob("*.EBD", ebd_dir);
ebd_path = [ebd_path1; ebd_path2];

sbd_path1 = Glob.glob("*.sbd", sbd_dir);
sbd_path2 = Glob.glob("*.SBD", sbd_dir);
sbd_path = [sbd_path1; sbd_path2];

tbd_path1 = Glob.glob("*.tbd", tbd_dir);
tbd_path2 = Glob.glob("*.TBD", tbd_dir);
tbd_path = [tbd_path1; tbd_path2];

# number of data files
ndatafiles = length(dbd_path);

# extracting list of cache files (using comprehension)
cac_list = [cac_path[i][end-11:end-4] for i in 1:length(cac_path)];

# extract the header of the dbd files
ndbdheaderlines = 14;
dbdheader = Array{String}(undef,ndatafiles,ndbdheaderlines);
ndbdsensors = Array{Int}(undef,ndatafiles);
dbdcacname = Array{String}(undef,ndatafiles);
cacind = Array{Int}(undef,ndatafiles);

# extract the header array for all the data files (dbdheader), the number of sensors for each data file (ndbdsensors), and the CAC name associated with each data file (dbdcacname)
for i = 1:ndatafiles
    #display(ndbdsensors)
    dbdio = open(dbd_path[i], "r")
    #display(dbd_path[i])
    for j = 1:ndbdheaderlines
        dbdheader[i,j] = readline(dbdio);
    end
    close(dbdio)

    # determine the number of DBD sensors & the cache file name for the DBD file
    ndbdsensors[i] = parse(Int64, dbdheader[i,10][end-4:end]);
    dbdcacname[i] = dbdheader[i,13][end-7:end];

    # calculate the index of dbdcacname in the list of cache files
    cacind[i] = findall(dbdcacname[i] .== cac_list)[1];
end

# construct sensorlist
sensorlist = SensorList[];
for i = 1:ndatafiles
    display(i)
    # load each line in the cache file
    cacdata = Array{String}(undef,ndbdsensors[i]);
    cacarray = Array{String}(undef,ndbdsensors[i],7);
    open(cac_path[cacind][i], "r") do cacio
        cacdata = readlines(cacio, keep = false);
    end

    # split the strings into arrays of values for each sensor
    cacdataarr = split.(cacdata," ", keepempty=false);

    # construct 2-D string array of dbd header
    for j = 1:length(cacdataarr)
        cacarray[j,:] = cacdataarr[j];
    end

    # change T to true and F to false in transmitted boolean array
    cacarray[:,2] = replace(x -> x=="T" ? "true" : "false", cacarray[:,2]);
    transmitted = parse.(Bool, cacarray[:,2]);

    # sensor number
    sensnum = parse.(Int, cacarray[:,3]);

    # sensor index number
    indexnum = parse.(Int, cacarray[:,4]);

    # number of bytes per entry
    nbytes = parse.(Int, cacarray[:,5]);

    name = cacarray[:,6];
    unit = cacarray[:,7];

    push!(sensorlist,SensorList(transmitted, sensnum, indexnum, nbytes, name, unit));
end 

#mutable struct data_struct
#    eval(Symbol(sensorlist[1].name[1], "::Array{Float64}"))
#end

