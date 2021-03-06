"""
    transform(::BoxCox, 饾惐)

Transform an array using Box-Cox method.  The power parameter `位` is derived
from maximizing a log-likelihood estimator.
If the array contains any non-positive values then a `DomainError` is thrown.
This can be avoided by providing the shift parameter `伪` to make all values
positive.

# Keyword arguments
- `伪`: added to all values in 饾惐 before transformation. Default = 0.
- `scaled`: scale transformation results.  Default = false.
"""
function transform(m::BoxCox, 饾惐; kwargs...)
    位, details = lambda(m, 饾惐; kwargs...)
    @debug "BoxCox: estimated lambda = $位"
    return transform(m, 饾惐, 位; kwargs...)
end

"""
    transform(::BoxCox, 饾惐, 位; 伪 = 0)

Transform an array using Box-Cox method with the provided power parameter `位`.
If the array contains any non-positive values then a `DomainError` is thrown.

# Keyword arguments
- `伪`: added to all values in 饾惐 before transformation. Default = 0.
- `scaled`: scale transformation results.  Default = false.
"""
function transform(::BoxCox, 饾惐, 位; 伪 = 0, scaled = false, kwargs...)
    if 伪 != 0
        饾惐 .+= 伪
    end
    any(饾惐 .<= 0) && throw(DomainError(
        "Data must be positive and ideally greater than 1. " *
        "You may specify 伪 argument(shift). "))
    if scaled
        gm = geomean(饾惐)
        return @. 位 鈮? 0 ? gm * log(饾惐) : (饾惐 ^ 位 - 1) / (位 * gm ^ (位 - 1))
    else
        return @. 位 鈮? 0 ? log(饾惐) : (饾惐 ^ 位 - 1) / 位
    end
end

"""
    lambda(::BoxCox, 饾惐; interval = (-2.0, 2.0), method = :geomean)

Calculate lambda from an array using a log-likelihood estimator.

# Keyword arguments
- `method`: either :geomean or :normal
- you can also pass any other keyword arguments accepted by `Optim.optimize` function e.g. `abs_tol`

See also: [`log_likelihood`](@ref)
"""
function lambda(m::BoxCox, 饾惐; interval = (-2.0, 2.0), kwargs...)
    i1, i2 = interval
    res = optimize(位 -> -log_likelihood(m, 饾惐, 位; kwargs...), i1, i2)
    return (value=minimizer(res), details=res)
end

"""
    log_likelihood(::BoxCox, 饾惐, 位; method = :geomean)

Return log-likelihood for the given array and lambda.

For method `:geomean`:
```
    -N / 2.0 * log(2 * 蟺 * 蟽虏 / gm ^ (2 * (位 - 1)) + 1)
```

For method `:normal`:
```
    -N / 2.0 * log(蟽虏) + (位 - 1) * sum(log.(饾惐))
```
"""
function log_likelihood(m::BoxCox, 饾惐, 位; method = :geomean, kwargs...)
    N = length(饾惐)
    饾惒 = transform(m, float.(饾惐), 位)
    蟽虏 = var(饾惒, corrected = false)
    gm = geomean(饾惐)
    if method == :geomean
        -N / 2.0 * log(2 * 蟺 * 蟽虏 / gm ^ (2 * (位 - 1)) + 1)
    elseif method == :normal
        -N / 2.0 * log(蟽虏) + (位 - 1) * sum(log.(饾惐))
    else
        throw(ArgumentError("Incorrect method. Please specify :geomean or :normal."))
    end
end
