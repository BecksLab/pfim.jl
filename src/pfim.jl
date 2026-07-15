module pfim

# Dependencies
using DataFrames
using SpeciesInteractionNetworks

include(joinpath("lib", "downsample.jl"))
export downsample_network

include(joinpath("lib", "pfim.jl"))
export PFIM

end # module pfim
