
"""
    PFIM(trait_data::DataFrame, feeding_rules::DataFrame; kwargs...)

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
        - `:network` → `SpeciesInteractionNetwork`
        - `:matrix` → adjacency matrix (`Bool`)
        - `:edgelist` → vector of `(resource, consumer)` tuples

    - `downsample::Bool = false`:
        Whether to apply probabilistic downsampling to reduce link density.

    - `y::Float64 = 2.5`:
        Downsampling parameter controlling expected number of links per consumer.

    # Returns
    Depends on `return_type`:

    - `:network` → `SpeciesInteractionNetwork`
    - `:matrix` → `Matrix{Bool}` of size `S × S`
    - `:edgelist` → `Vector{Tuple{Any, Any}}` of `(resource, consumer)` pairs

    # Details
    Trait matching is performed across all specified trait types. For each taxon pair,
    the number of satisfied feeding rules is counted. Interactions are retained if this
    count meets or exceeds the `certainty_req` threshold.

    If provided, the numerical size rule is applied in addition to categorical matching.

    Downsampling is applied after interaction inference and probabilistically removes
    links based on consumer degree.
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

    # data checks
    if !(certainty_req == :all || certainty_req isa Int)
        error("certainty_req must be :all or Int")
    end

    # --- derive trait types from rules ---
    if isnothing(trait_types)
        trait_types = Symbol.(unique(feeding_rules.trait_type_resource))
    end

    S = nrow(trait_data)
    int_matrix = zeros(Bool, S, S)

    # --- certainty threshold ---
    threshold = certainty_req == :all ? length(trait_types) : certainty_req

    # --- build interaction matrix ---
    for cons in 1:S
        for res in 1:S

            if !allow_self && cons == res
                continue
            end

            consumer = trait_data[cons, :]
            resource = trait_data[res, :]

            tally = 0

            for trait in trait_types
                consumer_trait = consumer[trait]
                resource_trait = resource[trait]

                allowed_resources =
                    feeding_rules[
                        feeding_rules.trait_type_consumer .== String(trait) .&&
                        feeding_rules.trait_consumer .== consumer_trait,
                        :trait_resource
                    ]

                if resource_trait ∈ allowed_resources
                    tally += 1
                end
            end

            # --- numeric size rule ---
            if !isnothing(size_col) && !isnothing(num_size_rule)
                res_size = resource[size_col]
                con_size = consumer[size_col]

                if num_size_rule(res_size, con_size) == 0
                    continue
                end
            end

            if tally >= threshold
                int_matrix[cons, res] = 1
            end
        end
    end

    taxa = Symbol.(trait_data[:, taxon_col])

    # --- downsampling step ---
    # essentially updates int_matrix with downsampled (Binary) version
    if downsample
        int_matrix = _downsample(int_matrix, taxa, y)
    end

    # --- output handling ---
    if return_type == :matrix
        return int_matrix

    elseif return_type == :edgelist
        edges = [
            (taxa[res], taxa[cons])
            for cons in 1:S, res in 1:S
            if int_matrix[cons, res] == 1
        ]
        return edges

    elseif return_type == :network
        nodes = Unipartite(Symbol.(taxa))
        edges = Binary(int_matrix)
        return SpeciesInteractionNetwork(nodes, edges)

    else
        error("return_type must be :matrix, :edgelist, or :network")
    end
end