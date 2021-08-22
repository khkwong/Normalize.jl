module Normalize

using Optim: optimize, minimizer
using Statistics: mean, var
using StatsBase: geomean, Histogram
using Distributions: fit, Normal, pdf

include("types.jl")
include("chisq.jl")
include("boxcox.jl")
include("yeojohnson.jl")

end
