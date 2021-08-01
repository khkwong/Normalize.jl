module Normalize

using Optim: optimize, minimizer
using Statistics: mean, var
using StatsBase: geomean
using Distributions

include("types.jl")
include("boxcox.jl")
include("chisq.jl")

end
