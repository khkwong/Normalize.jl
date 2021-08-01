module Normalize

using Optim: optimize, minimizer
using Statistics: mean, var
using StatsBase: geomean, Histogram
using Distributions: fit, Normal, pdf

include("types.jl")
include("boxcox.jl")
include("chisq.jl")

end
