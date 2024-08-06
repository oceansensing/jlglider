using DataFrames, CSV

jm22df = DataFrame("Time" => jm22.gpst, "Longtitude" => jm22.gpslon, "Latitude" => jm22.gpslat);
lbe22df = DataFrame("Time" => lbe22.gpst, "Longtitude" => lbe22.gpslon, "Latitude" => lbe22.gpslat);
jm23df = DataFrame("Time" => jm23.gpst, "Longtitude" => jm23.gpslon, "Latitude" => jm23.gpslat);

CSV.write("SEA064_JM2022_lonlat.csv", jm22df);
CSV.write("SEA064_LBE2022_lonlat.csv", lbe22df);
CSV.write("SEA064_JM2023_lonlat.csv", jm23df);
