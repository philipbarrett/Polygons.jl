using Polygons
using Base.Test

tests = ["definitions", "addition" ]
for t in tests
    include("$(t).jl")
end
