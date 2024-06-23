workdir = "/Users/gong/GitHub/jlglider/seaexplorer"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir);
end

include("write_glider_data.jl")
write_glider_data(glider = "SEA064", mission = 48, csvdir = "./");