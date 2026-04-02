
"""
   PFIM(trait_data::DataFrame, feeding_rules::DataFrame; y::Float64 = 2.5, downsample::Bool = true)

    Takes a data frame and implements the feeding rules to determine the
    feasibility of links between species. As well as applying the link
    distribution downsampling approach.
    
    #### References
    
    Shaw, Jack O., Alexander M. Dunhill, Andrew P. Beckerman, Jennifer A.
    Dunne, and Pincelli M. Hull. 2024. “A Framework for Reconstructing 
    Ancient Food Webs Using Functional Trait Data.” 
    https://doi.org/10.1101/2024.01.30.578036.
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

    taxa = trait_data[:, taxon_col]

    # --- downsampling step ---
    if downsample
        int_matrix = _downsample(int_matrix, y)
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