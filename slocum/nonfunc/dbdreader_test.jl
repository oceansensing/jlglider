using PyCall, Dates, NaNMath
dbdreader = pyimport("dbdreader");

function pyrol2jlcol(invar::Matrix{Float64})
    return reverse(rotr90(invar), dims = 2);
end

datadir = "/Users/gong/Research/electa-20221103-passengers";
t0 = DateTime("2022-10-01");
tN = DateTime("2022-12-01");

#dbdSylvia = dbdreader.DBD("/Users/gong/GitHub/sylvia-20180501-maracoos/data/DBD/01740010.DBD", cacheDir="/Users/gong/GitHub/sylvia-20180501-maracoos/data/cache/");
#dbdElecta = dbdreader.DBD(datadir * "/DBD/00970000.dbd", cacheDir = datadir * "/cache/");
#ebdElecta = dbdreader.DBD(datadir * "/EBD/00970000.ebd", cacheDir = datadir * "/cache/");

#dbdElecta = dbdreader.MultiDBD(datadir * "/DBD/0*.dbd", cacheDir = datadir * "/cache/");
#ebdElecta = dbdreader.MultiDBD(datadir * "/EBD/0*.ebd", cacheDir = datadir * "/cache/");

debdElecta = dbdreader.MultiDBD(pattern = datadir * "/[DE]BD/0*.[de]bd", complement_files = true, cacheDir = datadir * "/cache/");
dbdvars = debdElecta.parameterNames["eng"];
ebdvars = debdElecta.parameterNames["sci"];

#tctd, cond, temp, pres, m_de_oil_vol = debdElecta.get_CTD_sync("m_de_oil_vol");

sci_m_present_time = pyrol2jlcol(debdElecta.get("sci_m_present_time"));
sci_water_temp = pyrol2jlcol(debdElecta.get("sci_water_temp"));
sci_water_cond = pyrol2jlcol(debdElecta.get("sci_water_cond"));
sci_water_pressure = pyrol2jlcol(debdElecta.get("sci_water_pressure")); 

trange = datetime2unix.([t0; tN]);

tall = unique(union(sci_water_temp[:,1], sci_water_cond[:,1], sci_water_pressure[:,1]));
#tctd = unique(sci_m_present_time[1:end,1]);
tind = findall(trange[1] .< tall .< trange[end]);
tctd = tall[tind];

t1 = unix2datetime(floor(NaNMath.minimum(tctd)));
t2 = unix2datetime(ceil(NaNMath.maximum(tctd)));
tt = collect(t1:Second(1):t2);
