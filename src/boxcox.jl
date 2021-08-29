"""
    transform(::BoxCox, 𝐱)

Transform an array using Box-Cox method.  The power parameter `λ` is derived
from maximizing a log-likelihood estimator.
If the array contains any non-positive values then a `DomainError` is thrown.
This can be avoided by providing the shift parameter `α` to make all values
positive.

# Keyword arguments
- `α`: added to all values in 𝐱 before transformation. Default = 0.
- `scaled`: scale transformation results.  Default = false.
"""
function transform(m::BoxCox, 𝐱; kwargs...)
    λ, details = lambda(m, 𝐱; kwargs...)
    @debug "BoxCox: estimated lambda = $λ"
    return transform(m, 𝐱, λ; kwargs...)
end

"""
    transform(::BoxCox, 𝐱, λ; α = 0)

Transform an array using Box-Cox method with the provided power parameter `λ`.
If the array contains any non-positive values then a `DomainError` is thrown.

# Keyword arguments
- `α`: added to all values in 𝐱 before transformation. Default = 0.
- `scaled`: scale transformation results.  Default = false.
"""
function transform(::BoxCox, 𝐱, λ; α = 0, scaled = false, kwargs...)
    if α != 0
        𝐱 .+= α
    end
    any(𝐱 .<= 0) && throw(DomainError(
        "Data must be positive and ideally greater than 1. " *
        "You may specify α argument(shift). "))
    if scaled
        gm = geomean(𝐱)
        return @. λ ≈ 0 ? gm * log(𝐱) : (𝐱 ^ λ - 1) / (λ * gm ^ (λ - 1))
    else
        return @. λ ≈ 0 ? log(𝐱) : (𝐱 ^ λ - 1) / λ
    end
end

"""
    lambda(::BoxCox, 𝐱; interval = (-2.0, 2.0), method = :geomean)

Calculate lambda from an array using a log-likelihood estimator.

# Keyword arguments
- `method`: either :geomean or :normal
- you can also pass any other keyword arguments accepted by `Optim.optimize` function e.g. `abs_tol`

See also: [`log_likelihood`](@ref)
"""
function lambda(m::BoxCox, 𝐱; interval = (-2.0, 2.0), kwargs...)
    i1, i2 = interval
    res = optimize(λ -> -log_likelihood(m, 𝐱, λ; kwargs...), i1, i2)
    return (value=minimizer(res), details=res)
end

"""
    log_likelihood(::BoxCox, 𝐱, λ; method = :geomean)

Return log-likelihood for the given array and lambda.

For method `:geomean`:
```
    -N / 2.0 * log(2 * π * σ² / gm ^ (2 * (λ - 1)) + 1)
```

For method `:normal`:
```
    -N / 2.0 * log(σ²) + (λ - 1) * sum(log.(𝐱))
```
"""
function log_likelihood(m::BoxCox, 𝐱, λ; method = :geomean, kwargs...)
    N = length(𝐱)
    𝐲 = transform(m, float.(𝐱), λ)
    σ² = var(𝐲, corrected = false)
    gm = geomean(𝐱)
    if method == :geomean
        -N / 2.0 * log(2 * π * σ² / gm ^ (2 * (λ - 1)) + 1)
    elseif method == :normal
        -N / 2.0 * log(σ²) + (λ - 1) * sum(log.(𝐱))
    else
        throw(ArgumentError("Incorrect method. Please specify :geomean or :normal."))
    end
end
