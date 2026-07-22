module pfim

# Dependencies
using DataFrames
using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics

include(joinpath("lib", "downsample.jl"))
export downsample_network

include(joinpath("lib", "pfim.jl"))
export PFIM

end # module pfim
