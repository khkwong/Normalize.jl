using Normalize
using Test

@testset "Normalize.jl" begin
    include("test_boxcox.jl")
    include("test_yeojohnson.jl")
end
