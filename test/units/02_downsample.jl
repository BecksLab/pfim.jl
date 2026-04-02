using DataFrames
using pfim
using SpeciesInteractionNetworks
using Test

# -------------------------------
# 1. Downsampling works
# -------------------------------
@testset "Downsampling: basic structure" begin
    trait_data = DataFrame(
        species = [:a, :b, :c],
        trait1 = ["x", "y", "x"]
    )

    feeding_rules = DataFrame(
        trait_type_resource = ["trait1"],
        trait_resource = ["x"],
        trait_type_consumer = ["trait1"],
        trait_consumer = ["y"]
    )

    net = PFIM(trait_data, feeding_rules;
        return_type = :network,
        downsample = true
    )

    @test isa(net, SpeciesInteractionNetwork)
    @test richness(net) == 3
    @test species(net) == [:a, :b, :c]
end