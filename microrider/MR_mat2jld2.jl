#localpath = ENV["HOME"] * "/Users/gong/GitHub/jlglider/microrider"
#if (localpath in LOAD_PATH) == false
#    push!(LOAD_PATH, localpath);
#end

include("MR_io.jl")
using .MR_io
using Glob, MAT, JLD2

MR_io.MR_mat2jld2("NORSE", "JM", 2023)
