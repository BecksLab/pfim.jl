using Test

@testset "pfim" begin
    include("units/00_allgood.jl")
    include("units/01_metaweb.jl")
    include("units/02_downsample.jl")
end
