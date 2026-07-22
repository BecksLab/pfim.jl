"""
Downsampling.jl
-------------------------
Provides topological network downsampling routines based on species link 
distributions using the scaling framework from Roopnarine (2006).

### Matrix Convention:
- **Rows (i):** Consumers / Predators (who is doing the eating)
- **Columns (j):** Resources / Prey (who is being eaten)
- `matrix[i, j] = 1` indicates that consumer `i` eats resource `j`.

### Ecological Mechanics:
The algorithm scales a link's survival probability based on the consumer's total 
generality (out-degree, computed along rows). Highly generalist consumers retain links 
with higher log-odds baseline weights, while specialists or reduced webs prune down 
exponentially relative to a network scale constraint \$E\$.
"""

using LinearAlgebra
using Statistics

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

"""
    _get_downsample_metrics(mat::AbstractMatrix{Bool}, S::Int) -> Tuple{Int, Float64}

Calculate internal metrics used to monitor network degradation during downsampling.
"""
function _get_downsample_metrics(mat::AbstractMatrix{Bool}, S::Int)
    pred_has_links = vec(sum(mat, dims=2) .> 0) # Active consumer rows
    prey_has_links = vec(sum(mat, dims=1) .> 0) # Active resource columns
    active = sum(pred_has_links .|| prey_has_links)
    
    co = sum(mat) / (S^2)
    return active, co
end

"""
    _single_downsample_step(mat::AbstractMatrix{Bool}, y::Float64) -> Matrix{Bool}

Perform a single-step probabilistic link pruning on the network.
Implements the core scaling logic from Roopnarine (2006).
"""
function _single_downsample_step(mat::AbstractMatrix{Bool}, y::Float64)
    S = size(mat, 1)
    
    generality_vector = vec(sum(mat, dims=2)) 

    E = exp(log(S) * (y - 1) / y)
    link_dist = exp.(generality_vector ./ E)

    prob_matrix = zeros(Float64, S, S)
    for i in 1:S
        prey_indices = findall(mat[i, :])
        prob_matrix[i, prey_indices] .= link_dist[i]
    end

    maxval = maximum(prob_matrix)
    if maxval > 0 && isfinite(maxval)
        prob_matrix ./= maxval
    else
        prob_matrix .= 0.0
    end
    prob_matrix = clamp.(prob_matrix, 0.0, 1.0)

    # FIX: Mask with original matrix using bitwise AND to avoid generating fake links
    return mat .& (rand(S, S) .<= prob_matrix)
end

# Categorical sampler for weighted probabilities
function _rand_categorical(p::Vector{Float64})
    r = rand()
    cp = 0.0
    for i in 1:length(p)
        cp += p[i]
        if r <= cp
            return i
        end
    end
    return length(p)
end

# ==============================================================================
# Public API Function
# ==============================================================================

"""
    downsample_network(int_matrix::AbstractMatrix{Bool}, y::Float64; kwargs...) -> Matrix{Bool}

Downsample a food web's interaction matrix based on species link distributions.

Supports both standard single-step probabilistic pruning and iterative pruning 
targeted to match a specific network connectance (Co).
"""
function downsample_network(
    int_matrix::AbstractMatrix{Bool}, 
    y::Float64;
    target_co::Union{Nothing, Float64} = nothing,
    min_spp_prop::Float64 = 0.5,
    max_iter::Int = 50
)
    S = size(int_matrix, 1)

    # --- Case 1: Single-Step Downsampling ---
    if isnothing(target_co)
        return _single_downsample_step(int_matrix, y)
    end

    # --- Case 2: Iterative Connectance-Targeted Downsampling ---
    current_matrix = copy(int_matrix)
    min_species = ceil(Int, S * min_spp_prop)
    _, current_co = _get_downsample_metrics(current_matrix, S)

    if current_co <= target_co
        @warn "Initial connectance ($current_co) is already <= target ($target_co). Returning original network."
        return current_matrix
    end

    # Track the closest network encountered
    best_matrix = copy(current_matrix)
    best_co = current_co
    best_diff = abs(current_co - target_co)

    iter = 0

    while current_co > target_co && iter < max_iter
        iter += 1

        links = findall(current_matrix)
        if isempty(links)
            break
        end

        # Calculate Roopnarine (2006) weights for existing links
        generality_vector = vec(sum(current_matrix, dims=2))
        E = exp(log(S) * (y - 1) / y)
        link_dist = exp.(generality_vector ./ E)

        # Compute prune weights
        link_prune_weights = Vector{Float64}(undef, length(links))
        for (idx, link) in enumerate(links)
            i = link[1]
            p_retain = link_dist[i]
            link_prune_weights[idx] = 1.0 - p_retain + 1e-6
        end

        prob_dist = link_prune_weights ./ sum(link_prune_weights)

        # Select one link to prune
        chosen_idx = _rand_categorical(prob_dist)
        target_link = links[chosen_idx]

        temp_matrix = copy(current_matrix)
        temp_matrix[target_link] = false

        active_spp, next_co = _get_downsample_metrics(temp_matrix, S)

        # Reject removals that violate species retention
        if active_spp < min_species
            continue
        end

        # Accept the removal
        current_matrix = temp_matrix
        current_co = next_co

        # Update best network if this one is closer to the target
        current_diff = abs(current_co - target_co)

        if current_diff < best_diff
            best_diff = current_diff
            best_co = current_co
            best_matrix = copy(current_matrix)
        end
    end

    if iter == max_iter && current_co > target_co
        @warn "Reached max iterations ($max_iter) before dropping below target connectance. Closest connectance found: $best_co"
    end

    return best_matrix
end