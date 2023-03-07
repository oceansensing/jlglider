module slocum_types

using Dates

mutable struct SensorListTest2
    Symbol("transmitted::Array{Bool}")
    sensnum::Array{UInt16}
    indexnum::Array{Int16}
    nbytes::Array{UInt8}
    name::Array{String}
    unit::Array{String}
end

mutable struct GliderDataStruct
    for ii = 1:length(sensorlist[1].name)
        #sensorlist[1].name[ii]
        Symbol(sensorlist[1].name[2], "::Array{Float64}")
    end
end

mutable struct TestType
    Symbol(sensorlist[1].name[1], "::Float64")
end


end