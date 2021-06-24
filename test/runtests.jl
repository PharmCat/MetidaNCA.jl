using Test
using DataFrames, CSV

path     = dirname(@__FILE__)
io       = IOBuffer();
pkdata2  = CSV.File(joinpath(path, "csv", "pkdata2.csv")) |> DataFrame

@testset "MetidaNCA.jl" begin

ds = MetidaNCA.pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation])

sbj = MetidaNCA.nca!(ds[1])
end
