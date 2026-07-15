using DataFrames
using pfim
using SpeciesInteractionNetworks
using Test

# -------------------------------
# 1. Downsampling works
# -------------------------------
@testset "Downsampling: basic structure" begin
    # 1. Define species with distinct traits to create asymmetric diets
    trait_data = DataFrame(
        species = [:a, :b, :c, :d, :e],
        trait1  = ["x", "y", "z", "w", "x"]  # b is 'y', d is 'w'
    )

    # 2. Feeding rules that divide generalists and specialists:
    # - Consumer 'y' (:b) eats 'x' (:a, :e) and 'z' (:c)  --> Generality = 3
    # - Consumer 'w' (:d) eats ONLY 'x' (:a, :e)          --> Generality = 2
    feeding_rules = DataFrame(
        trait_type_resource = ["trait1", "trait1", "trait1"],
        trait_resource      = ["x", "z", "x"],
        trait_type_consumer = ["trait1", "trait1", "trait1"],
        trait_consumer      = ["y", "y", "w"]
    )

    # --- Step 1: Run PFIM WITHOUT downsampling ---
    net_full = PFIM(trait_data, feeding_rules;
        return_type = :network,
        downsample = false
    )
    full_matrix = net_full.edges.edges
    initial_link_count = sum(full_matrix) # Should be 5 links

    # --- Step 2: Run PFIM WITH downsampling ---
    #Random.seed!(66)
    net_pruned = PFIM(trait_data, feeding_rules;
        return_type = :network,
        downsample = true,
        y = 1.1  # Aggressive scaling to ensure lower-generality links drop
    )

    pruned_matrix = net_pruned.edges.edges
    pruned_link_count = sum(pruned_matrix)

    # --- Step 3: Assertions ---
    @test isa(net_pruned, SpeciesInteractionNetwork)
    @test richness(net_pruned) == 5
    @test species(net_pruned) == [:a, :b, :c, :d, :e]

    # These will now pass!
    @test pruned_link_count < initial_link_count
    @test all(pruned_matrix .<= full_matrix)
end

# -------------------------------
# 2. External Downsampling & Target Connectance
# -------------------------------
@testset "Downsampling: target connectance & external use" begin
    # Create a dense 5x5 matrix (connectance = 0.8)
    # S = 5, S^2 = 25 potential links. 20 active links.
    taxa = [:sp1, :sp2, :sp3, :sp4, :sp5]
    dense_matrix = Bool[
        1 1 1 1 0;
        1 1 1 1 1;
        0 1 1 1 1;
        1 1 0 1 1;
        1 1 1 0 1
    ]
    
    # Verify initial connectance is 20/25 = 0.8
    init_co = sum(dense_matrix) / 25
    @test init_co == 0.84

    # A: Test basic external single-step downsample (no target_co)
    single_down = downsample_network(dense_matrix, taxa, 2.5)
    @test size(single_down) == (5, 5)
    @test eltype(single_down) == Bool

    # B: Test targeting a specific lower connectance (e.g., 0.4)
    target_co = 0.4
    pruned_matrix = downsample_network(dense_matrix, taxa, 2.5; target_co = target_co)
    
    final_co = sum(pruned_matrix) / 25
    @test final_co <= target_co
    
    # C: Test Defensive Guardrails (extremely low target connectance)
    # If we target 0.05 (which requires dropping almost all links), 
    # we should halt before we violate min_spp_prop (default 0.5, meaning 3 active species)
    guarded_matrix = downsample_network(
        dense_matrix, 
        taxa, 
        2.5; 
        target_co = 0.04, 
        min_spp_prop = 0.6  # Needs at least 3 species to keep a link
    )
    
    # Count species that still have interactions
    active_spp = sum(sum(guarded_matrix, dims=1) .> 0 .|| sum(guarded_matrix, dims=2)' .> 0)
    
    # We should have respected the guardrail of retaining at least 3 active species
    @test active_spp >= 3
end