module SPTestMetaweb

using CSV
using DataFrames
using pfim
using SpeciesInteractionNetworks
using Test

# -------------------------------
# Load test data
# -------------------------------
feeding_rules = DataFrame(CSV.File("units/data/feeding_rules.csv"))
traits = DataFrame(CSV.File("units/data/trait.csv"))
traits = DataFrame([col => string.(traits[!, col]) for col in names(traits)])

known_list = DataFrame(CSV.File("units/data/interactions.csv"))
# Keep known_edges structured natively as (Resource, Consumer) from the file
known_edges = [
    (Symbol(r.resource), Symbol(r.consumer)) for r in eachrow(known_list)
]

# Create a mapping of known edges to our new corrected convention: (Consumer, Resource)
expected_consumer_resource_edges = [
    (c, r) for (r, c) in known_edges
]

# -------------------------------
# 1. Runs without error
# -------------------------------
@testset "PFIM runs without error" begin
    @test_nowarn pfim.PFIM(traits, feeding_rules; downsample=false)
end

# -------------------------------
# 2. Network output matches expected
# -------------------------------
@testset "Network output matches expected interactions" begin
    net = pfim.PFIM(traits, feeding_rules; return_type=:network, downsample=false)
    pfim_int = interactions(net)

    # SpeciesInteractionNetworks uses (from, to, value) -> (Consumer, Resource, true)
    @test sort(pfim_int) == sort([(c, r, true) for (c, r) in expected_consumer_resource_edges])
end

# -------------------------------
# 3. Edgelist output matches expected
# -------------------------------
@testset "Edgelist output matches expected interactions" begin
    edges = pfim.PFIM(traits, feeding_rules; return_type=:edgelist, downsample=false)

    # FIX: Now matching against our corrected (Consumer, Resource) layout
    @test sort(edges) == sort(expected_consumer_resource_edges)
end

# -------------------------------
# 4. Matrix output is consistent
# -------------------------------
@testset "Matrix output is consistent with edgelist" begin
    mat = pfim.PFIM(traits, feeding_rules; return_type=:matrix, downsample=false)
    edges = pfim.PFIM(traits, feeding_rules; return_type=:edgelist, downsample=false)

    taxa = Symbol.(traits.species)
    edge_set = Set(edges)

    # FIX: Adjusted reconstructed tuple to (Consumer, Resource) to mirror 
    # the underlying matrix rows (cons) and columns (res).
    reconstructed = Set([
        (taxa[cons], taxa[res])
        for cons in axes(mat,1), res in axes(mat,2)
        if mat[cons, res]
    ])

    @test edge_set == reconstructed
end

# -------------------------------
# 5. Certainty requirement behaviour
# -------------------------------
@testset "certainty_req filters interactions correctly" begin
    strict = pfim.PFIM(traits, feeding_rules;
        certainty_req=:all,
        return_type=:edgelist,
        downsample=false
    )

    relaxed = pfim.PFIM(traits, feeding_rules;
        certainty_req=1,
        return_type=:edgelist,
        downsample=false
    )

    @test length(strict) <= length(relaxed)
end

# -------------------------------
# 6. Missing taxon column throws error
# -------------------------------
@testset "Missing taxon column throws error" begin
    bad_traits = deepcopy(traits)
    rename!(bad_traits, :species => :wrong)

    @test_throws Exception pfim.PFIM(
        bad_traits,
        feeding_rules;
        taxon_col=:species,
        downsample=false
    )
end

# -------------------------------
# 7. Invalid certainty_req throws error
# -------------------------------
@testset "Invalid certainty_req errors" begin
    @test_throws Exception pfim.PFIM(
        traits,
        feeding_rules;
        certainty_req="invalid",
        downsample=false
    )
end

# -------------------------------
# 8. Numeric size rule behaviour
# -------------------------------
@testset "Numeric size rule works correctly" begin

    numeric_traits = DataFrame(
        species = Symbol.(["Lion", "Zebra", "Grass"]),
        diet = ["carnivore", "herbivore", "producer"],
        body_mass = [190.0, 300.0, 1.0]
    )

    simple_rules = DataFrame(
        trait_type_resource = ["diet"],
        trait_resource = ["herbivore"],
        trait_type_consumer = ["diet"],
        trait_consumer = ["carnivore"]
    )

    mass_rule = (res, con) -> con >= 0.5 * res ? 1 : 0

    edges = pfim.PFIM(
        numeric_traits,
        simple_rules;
        taxon_col=:species,
        size_col=:body_mass,
        num_size_rule=mass_rule,
        certainty_req=:all,
        return_type=:edgelist,
        downsample=false
    )

    # to (:Lion, :Zebra) [Consumer, Resource]
    @test Symbol.(("Lion", "Zebra")) in edges

    # Stricter rule
    strict_rule = (res, con) -> con > res ? 1 : 0

    edges_strict = pfim.PFIM(
        numeric_traits,
        simple_rules;
        taxon_col=:species,
        size_col=:body_mass,
        num_size_rule=strict_rule,
        certainty_req=:all,
        return_type=:edgelist,
        downsample=false
    )

    @test isempty(edges_strict)
end

end