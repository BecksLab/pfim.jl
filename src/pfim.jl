module pfim

# Dependencies
using DataFrames
using SpeciesInteractionNetworks

include(joinpath("lib", "downsample.jl"))

include(joinpath("lib", "pfim.jl"))
export PFIM

end # module pfim
