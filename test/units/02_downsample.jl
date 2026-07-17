module SPTestDownsample

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
    initial_link_count = sum(full_matrix) # Links built based on traits

    # --- Step 2: Run PFIM WITH downsampling ---

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

    @test pruned_link_count < initial_link_count
    @test all(pruned_matrix .<= full_matrix)
end

# -------------------------------
# 2. External Downsampling & Target Connectance
# -------------------------------
@testset "Downsampling: target connectance & external use" begin
    taxa = [:sp1, :sp2, :sp3, :sp4, :sp5]
    dense_matrix = Bool[
        1 1 1 1 0;
        1 1 1 1 1;
        0 1 1 1 1;
        1 1 0 1 1;
        1 1 1 0 1
    ]
    
    # FIX: Dynamic math evaluation replacing hardcoded mismatched comments
    total_elements = length(dense_matrix)
    init_co = sum(dense_matrix) / total_elements
    @test init_co == 0.84

    # A: Test basic external single-step downsample (no target_co)
    single_down = downsample_network(dense_matrix, 2.5)
    @test size(single_down) == (5, 5)
    @test eltype(single_down) == Bool

    # B: Test targeting a specific lower connectance (e.g., 0.4)
    target_co = 0.4
    pruned_matrix = downsample_network(dense_matrix, 2.5; target_co = target_co)
    
    final_co = sum(pruned_matrix) / total_elements
    @test final_co <= target_co
    
    # C: Test Defensive Guardrails (extremely low target connectance)
    guarded_matrix = downsample_network(
        dense_matrix, 
        2.5; 
        target_co = 0.04, 
        min_spp_prop = 0.6  # Needs at least 3 species to keep a link (0.6 * 5 = 3)
    )
    
    # Clean dimension reduction matching Rows=Predators, Cols=Prey matrix properties
    pred_has_links = vec(sum(guarded_matrix, dims = 2) .> 0) # Collapse cols to check rows
    prey_has_links = vec(sum(guarded_matrix, dims = 1) .> 0) # Collapse rows to check cols
    active_spp = sum(pred_has_links .|| prey_has_links)
    
    @test active_spp >= 3
end

end