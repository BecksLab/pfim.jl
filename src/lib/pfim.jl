"""
PFIM.jl
-----------------
Generates ancient and modern unipartite food webs using the Paleo Food Web 
Inference Model (PFIM).

### Model Background:
The PFIM is a trait-based trophic inference framework designed to reconstruct 
food webs, particularly for paleoecological or data-scarce communities where direct 
dietary observations are unavailable (e.g., fossil records). The model maps 
taxonomic traits to categorical feeding rules to establish "feasible" feeding 
interactions, which can then be paired with size constraints or downsampled 
to approximate "realized" interaction densities.

### Matrix Convention:
- **Rows (i):** Consumers / Predators (who is doing the eating)
- **Columns (j):** Resources / Prey (who is being eaten)
- `matrix[i, j] = 1` indicates that consumer `i` eats resource `j`.
- Fully aligned with SpeciesInteractionNetworks.jl and EcologicalNetworksDynamics.jl.
"""

"""
    PFIM(trait_data::DataFrame, feeding_rules::DataFrame; kwargs...) -> SpeciesInteractionNetwork / Matrix / Vector

Infer a trophic interaction network using trait-based feeding rules following
the Paleo Food Web Inference Model (PFIM).

This function evaluates all possible consumer–resource pairs and determines
interaction feasibility based on categorical trait matching rules and,
optionally, a numerical size constraint. Interactions are retained if they
meet a specified certainty threshold. An optional downsampling step can be
applied to match expected link distributions.

# Arguments
- `trait_data::DataFrame`:
    A data frame where each row represents a taxon and columns contain trait values.
- `feeding_rules::DataFrame`:
    A data frame defining allowed consumer–resource trait combinations. Must contain:
    - `trait_type_resource`
    - `trait_resource`
    - `trait_type_consumer`
    - `trait_consumer`

# Keyword Arguments
- `taxon_col::Symbol = :species`:
    Column name in `trait_data` containing taxon identifiers.
- `trait_types::Union{Nothing, Vector{Symbol}} = nothing`:
    Subset of trait columns to use. If `nothing`, trait types are inferred from `feeding_rules`.
- `size_col::Union{Nothing, Symbol} = nothing`:
    Column containing numerical size values for taxa.
- `num_size_rule::Union{Function, Nothing} = nothing`:
    Function defining predator–prey size feasibility. Must take `(resource_size, consumer_size)`
    and return `1` (feasible) or `0` (infeasible).
- `certainty_req::Union{Symbol, Int} = :all`:
    Number of trait rules required for an interaction:
    - `:all` → all trait types must match
    - `Int` → minimum number of matching trait types
- `allow_self::Bool = true`:
    Whether to allow self-interactions (cannibalism).
- `return_type::Symbol = :network`:
    Output format:
    - `:network` → `SpeciesInteractionNetwork` (from SpeciesInteractionNetworks.jl)
    - `:matrix` → adjacency matrix (`Bool`)
    - `:edgelist` → vector of `(consumer, resource)` tuples [Aligned to (Predator, Prey)]
- `downsample::Bool = false`:
    Whether to apply probabilistic downsampling to reduce link density.
- `y::Float64 = 2.5`:
    Downsampling parameter controlling expected number of links per consumer.

# Returns
- `return_type = :matrix`  → `Matrix{Bool}` of size `S × S`
- `return_type = :edgelist` → `Vector{Tuple{Symbol, Symbol}}` of `(predator, prey)` pairs
- `return_type = :network`  → `SpeciesInteractionNetwork` with `Unipartite` nodes and `Binary` edges.
"""
function PFIM(
    trait_data::DataFrame,
    feeding_rules::DataFrame;
    taxon_col::Symbol = :species,
    trait_types::Union{Nothing, Vector{Symbol}} = nothing,
    size_col::Union{Nothing, Symbol} = nothing,
    num_size_rule::Union{Function, Nothing} = nothing,
    certainty_req::Union{Symbol, Int} = :all,
    allow_self::Bool = true,
    return_type::Symbol = :network,
    downsample::Bool = false,
    y::Float64 = 2.5,
)
    # Data check
    if !(certainty_req == :all || certainty_req isa Int)
        error("certainty_req must be :all or Int")
    end

    if isnothing(trait_types)
        trait_types = Symbol.(unique(feeding_rules.trait_type_resource))
    end

    S = nrow(trait_data)
    int_matrix = zeros(Bool, S, S)
    threshold = certainty_req == :all ? length(trait_types) : certainty_req

    # Pre-compile feeding rules into a Set-lookup Dict
    rules_lookup = Dict{Tuple{Symbol, Any}, Set{Any}}()
    for r in eachrow(feeding_rules)
        # Standardise keys as Strings or Symbols to match data vectors safely
        key = (Symbol(r.trait_type_consumer), r.trait_consumer)
        if !haskey(rules_lookup, key)
            rules_lookup[key] = Set{Any}()
        end
        push!(rules_lookup[key], r.trait_resource)
    end

    # --- Build Interaction Matrix ---
    # Rows (cons) = Predators, Columns (res) = Prey
    for cons in 1:S
        for res in 1:S
            if !allow_self && cons == res
                continue
            end

            tally = 0

            for trait in trait_types
                consumer_trait = trait_data[cons, trait]
                resource_trait = trait_data[res, trait]

                # Fast O(1) dictionary key mapping
                key = (trait, consumer_trait)
                if haskey(rules_lookup, key) && (resource_trait ∈ rules_lookup[key])
                    tally += 1
                end
            end

            # Numeric size rule
            if !isnothing(size_col) && !isnothing(num_size_rule)
                if num_size_rule(trait_data[res, size_col], trait_data[cons, size_col]) == 0
                    continue
                end
            end

            if tally >= threshold
                int_matrix[cons, res] = true
            end
        end
    end

    taxa = Symbol.(trait_data[:, taxon_col])

    if downsample
        int_matrix = downsample_network(int_matrix, y)
    end

    # --- Output Handling ---
    if return_type == :matrix
        return int_matrix

    elseif return_type == :edgelist
        edges = [
            (taxa[cons], taxa[res]) 
            for cons in 1:S, res in 1:S
            if int_matrix[cons, res]
        ]
        return edges

    elseif return_type == :network
        nodes = Unipartite(taxa)
        edges = Binary(int_matrix)
        return SpeciesInteractionNetwork(nodes, edges)
    else
        error("return_type must be :matrix, :edgelist, or :network")
    end
end