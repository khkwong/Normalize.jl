"""
    transform(::YeoJohnson, 𝐱)

Transform an array using Yeo-Johnson method.  The power parameter λ is derived
from maximizing a log-likelihood estimator. 
"""
function transform(method::YeoJohnson, 𝐱; optim_args...)
    λ, details = lambda(method, 𝐱; optim_args...)
    @debug "YeoJohnson: estimated lambda = $λ"
    return transform(method, 𝐱, λ)
end

"""
    transform(::YeoJohnson, 𝐱, λ)

Transform an array using Yeo-Johnson method with the provided power parameter λ. 
"""
function transform(::YeoJohnson, 𝐱, λ) 
    𝐱′ = similar(𝐱, Float64)
    for (i, x) in enumerate(𝐱)
        if x >= 0
            𝐱′[i] = λ ≈ 0 ? log(x + 1) : ((x + 1)^λ - 1)/λ 
        else
            𝐱′[i] = λ ≈ 2 ? -log(-x + 1) : -((-x + 1)^(2 - λ) - 1) / (2 - λ)
        end
    end
    return 𝐱′
end

"""
    lambda(::YeoJohnson, 𝐱; interval = (-2.0, 2.0), optim_args...)

Calculate lambda from an array using a log-likelihood estimator.
Keyword arguments:
- interval: search interval
- optim_args: keyword arguments accepted by Optim.optimize function
See also: [`log_likelihood`](@ref)
"""
function lambda(method::YeoJohnson, 𝐱; interval = (-2.0, 2.0), optim_args...)
    i1, i2 = interval
    res = optimize(λ -> -log_likelihood(method, 𝐱, λ), i1, i2; optim_args...)
    return (value=minimizer(res), details=res)
end

"""
    log_likelihood(::YeoJohnson, 𝐱, λ)

Return log-likelihood for the given array and lambda.
"""
function log_likelihood(method::YeoJohnson, 𝐱, λ)
    N = length(𝐱)
    𝐲 = transform(method, float.(𝐱), λ)
    σ² = var(𝐲, corrected = false)
    c = sum(sign.(𝐱) .* log.(abs.(𝐱) .+ 1))
    llf = -N / 2.0 * log(σ²) + (λ - 1) * c
    @debug "YeoJohnson: λ = $λ => σ²=$σ², c=$c, llf=$llf"
    return llf
end
