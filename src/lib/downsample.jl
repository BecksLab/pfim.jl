# ==============================================================================
# Helper Functions
# ==============================================================================

# Helper to calculate active species and current connectance
"""
    _get_downsample_metrics(mat::AbstractMatrix{Bool}, S::Int) -> Tuple{Int, Float64}

Calculate internal metrics used to monitor network degradation during downsampling.

# Arguments
- `mat::AbstractMatrix{Bool}`: Binary adjacency matrix (size `S × S`).
- `S::Int`: Total number of species in the pool.

# Returns
- `active::Int`: The count of species that have at least one incoming or outgoing interaction.
- `co::Float64`: The network connectance calculated as:
  
  Co = L/S^2
  
  where L is the number of active links in `mat`.
"""
function _get_downsample_metrics(mat::AbstractMatrix{Bool}, S::Int)
    # Active species: has at least one incoming or outgoing link
    active = sum(sum(mat, dims=1) .> 0 .|| sum(mat, dims=2)' .> 0)
    co = sum(mat) / (S^2)
    return active, co
end

# Single-step probabilistic prune (Roopnarine aligned)
"""
    _single_downsample_step(mat::AbstractMatrix{Bool}, taxa::Vector, y::Float64) -> Matrix{Bool}

Perform a single-step probabilistic link pruning on the network.

This function implements the core probabilistic scaling logic from Roopnarine (2006). 
Link retention probability scale is parameterized by the generality (out-degree) of 
consumers and a scaling exponent y.

# Arguments
- `mat::AbstractMatrix{Bool}`: Binary adjacency matrix where rows are consumers and columns are resources.
- `taxa::Vector`: Species identifiers.
- `y::Float64`: Scaling exponent controlling the expected link distributions.

# Mathematical Details
For each consumer species i with generality r_i (number of resources consumed), the scaling benchmark E is:

E = \\exp\\left(\\frac{\\ln(S) ⋅ (y - 1)}{y}\\right)

The probability of retaining an interaction for consumer i is proportional to:

p_i = \\exp\\left(\\frac{r_i}{E}\\right)

These probabilities are projected onto the adjacency matrix, normalized by the maximum probability value, and finally subjected to independent Bernoulli trials via element-wise comparison with a random distribution.
"""
function _single_downsample_step(mat::AbstractMatrix{Bool}, taxa::Vector, y::Float64)
    S = length(taxa)
    
    # Generality of each consumer
    generality_vector = vec(sum(mat, dims=2)) 

    # Calculate link distributions based on Roopnarine (2006)
    E = exp(log(S) * (y - 1) / y)
    link_dist = exp.(generality_vector ./ E)

    # Populate probability matrix directly
    prob_matrix = zeros(Float64, S, S)
    for i in 1:S
        for j in 1:S
            if mat[i, j]
                prob_matrix[i, j] = link_dist[i]
            end
        end
    end

    # --- FIX: Normalize FIRST, then Clamp ---
    maxval = maximum(prob_matrix)
    if maxval > 0 && isfinite(maxval)
        prob_matrix ./= maxval
    else
        prob_matrix .= 0.0
    end
    prob_matrix = clamp.(prob_matrix, 0.0, 1.0) # Safe final step

    # Probabilistic draw
    random_draw_matrix = rand(S, S) .<= prob_matrix

    return random_draw_matrix
end

# ==============================================================================
# Public API Function
# ==============================================================================

"""
    downsample_network(int_matrix::AbstractMatrix{Bool}, taxa::Vector, y::Float64; kwargs...) -> Matrix{Bool}

Downsample a food web's interaction matrix based on species link distributions.

This function supports both standard single-step probabilistic pruning and iterative 
pruning targeted to match a specific network connectance (Co). 

# Arguments
- `int_matrix::AbstractMatrix{Bool}`: Binary adjacency matrix (size `S × S`).
- `taxa::Vector`: A list of unique species identifiers of length `S`.
- `y::Float64`: Structural scaling parameter (typically around 2.0 to 3.0) 
  controlling expected links per consumer.

# Keyword Arguments
- `target_co::Union{Nothing, Float64} = nothing`: 
  The desired target connectance (L/S^2). If `nothing` (default), the function runs a 
  single-step probabilistic pruning and exits. If a float is provided, it iteratively 
  pruning the network until this connectance threshold is met or exceeded.
- `min_spp_prop::Float64 = 0.5`: 
  Defensive safeguard. The minimum proportion of the initial species pool that must retain 
  at least one link (incoming or outgoing). If an iterative step would violate this limit, 
  downsampling is halted to protect network integrity.
- `max_iter::Int = 50`: 
  Defensive safeguard. Maximum number of iterations to run in targeted mode to prevent 
  infinite loops.

# Returns
- `Matrix{Bool}`: The downsampled binary adjacency matrix.

# Details & Safeguards
When `target_co` is supplied, the network is pruned iteratively because a single-pass 
Bernoulli draw cannot guarantee a precise global connectance target on discrete graphs.
During each iteration, the probabilities are dynamically updated based on the *updated* state of the network.

To prevent catastrophic web collapse (e.g., losing too many species or ending up with a 
completely empty matrix), two safety guardrails will trigger an early break:
1. **Species Conservation:** If the proportion of functionally active species drops below `min_spp_prop`, the process aborts and returns the matrix state from the previous iteration.
2. **Total Collapse:** If the total link count reaches zero, the process aborts.

# References
- Roopnarine, Peter D. 2006. “Extinction Cascades and Catastrophe in Ancient Food Webs.” 
  *Paleobiology* 32 (1): 1-19. https://www.jstor.org/stable/4096814
"""
function downsample_network(
    int_matrix::AbstractMatrix{Bool}, 
    taxa::Vector, 
    y::Float64;
    target_co::Union{Nothing, Float64} = nothing,
    min_spp_prop::Float64 = 0.5,
    max_iter::Int = 50
)
    S = length(taxa)
    
    # --- Case 1: Single-Step Downsampling ---
    if isnothing(target_co)
        return _single_downsample_step(int_matrix, taxa, y)
    end

    # --- Case 2: Iterative Connectance-Targeted Downsampling ---
    current_matrix = copy(int_matrix)
    min_species = ceil(Int, S * min_spp_prop)
    _, current_co = _get_downsample_metrics(current_matrix, S)

    if current_co <= target_co
        @warn "Initial connectance ($current_co) is already <= target ($target_co)."
        return current_matrix
    end

    iter = 0
    while current_co > target_co && iter < max_iter
        iter += 1
        next_matrix = _single_downsample_step(current_matrix, taxa, y)
        active_spp, next_co = _get_downsample_metrics(next_matrix, S)

        # Did we lose too many active species?
        if active_spp < min_species
            @warn "Downsampling halted (iter $iter): Species retention threshold violated ($active_spp/$S active, limit is $min_species)."
            break
        end

        # Did the network completely collapse?
        if sum(next_matrix) == 0
            @warn "Downsampling halted (iter $iter): Network collapsed to 0 links."
            break
        end

        current_matrix = next_matrix
        current_co = next_co
    end

    if iter == max_iter && current_co > target_co
        @warn "Reached max iterations ($max_iter) without hitting target connectance. Current Co: $current_co"
    end

    return current_matrix
end