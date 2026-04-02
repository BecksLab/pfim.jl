"""
    _downsample(network, matrix, y)

    Internal function to downsample a network based on species link 
    distributions.

    #### References
    Roopnarine, Peter D. 2006. “Extinction Cascades and Catastrophe in
    Ancient Food Webs.” Paleobiology 32 (1): 1-19. 
    https://www.jstor.org/stable/4096814.

"""
function _downsample(int_matrix, taxa, y::Float64)

    spp = Unipartite(taxa)
    edg = Binary(int_matrix)
    _N = SpeciesInteractionNetwork(spp, edg)

    S = length(taxa)

    link_dist = zeros(Float64, S)

    for i in eachindex(taxa)
        sp = taxa[i]
        r = generality(_N, sp)
        E = exp(log(S) * (y - 1) / y)
        link_dist[i] = exp(r / E)
    end

    prob_matrix = zeros(Float64, size(int_matrix))

    for i in axes(int_matrix, 1)
        for j in axes(int_matrix, 2)
            if int_matrix[i, j]
                prob_matrix[i, j] = link_dist[i]
            end
        end
    end

    prob_matrix ./= maximum(prob_matrix)

    nodes = Unipartite(taxa)
    edges = Probabilistic(prob_matrix)
    N_down = SpeciesInteractionNetwork(nodes, edges)
    N_down = randomdraws(N_down)

    return Matrix(N_down.edges.edges)
end
