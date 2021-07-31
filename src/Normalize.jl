module Normalize

using Optim: optimize, minimizer
using Statistics: mean, var
using StatsBase: geomean

include("types.jl")
include("boxcox.jl")

end
