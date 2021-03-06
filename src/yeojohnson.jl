"""
    transform(::YeoJohnson, 饾惐)

Transform an array using Yeo-Johnson method.  The power parameter 位 is derived
from maximizing a log-likelihood estimator. 
"""
function transform(method::YeoJohnson, 饾惐; optim_args...)
    位, details = lambda(method, 饾惐; optim_args...)
    @debug "YeoJohnson: estimated lambda = $位"
    return transform(method, 饾惐, 位)
end

"""
    transform(::YeoJohnson, 饾惐, 位)

Transform an array using Yeo-Johnson method with the provided power parameter 位. 
"""
function transform(::YeoJohnson, 饾惐, 位) 
    饾惐鈥? = similar(饾惐, Float64)
    for (i, x) in enumerate(饾惐)
        if x >= 0
            饾惐鈥瞇i] = 位 鈮? 0 ? log(x + 1) : ((x + 1)^位 - 1)/位 
        else
            饾惐鈥瞇i] = 位 鈮? 2 ? -log(-x + 1) : -((-x + 1)^(2 - 位) - 1) / (2 - 位)
        end
    end
    return 饾惐鈥?
end

"""
    lambda(::YeoJohnson, 饾惐; interval = (-2.0, 2.0), optim_args...)

Calculate lambda from an array using a log-likelihood estimator.
Keyword arguments:
- interval: search interval
- optim_args: keyword arguments accepted by Optim.optimize function
See also: [`log_likelihood`](@ref)
"""
function lambda(method::YeoJohnson, 饾惐; interval = (-2.0, 2.0), optim_args...)
    i1, i2 = interval
    res = optimize(位 -> -log_likelihood(method, 饾惐, 位), i1, i2; optim_args...)
    return (value=minimizer(res), details=res)
end

"""
    log_likelihood(::YeoJohnson, 饾惐, 位)

Return log-likelihood for the given array and lambda.
"""
function log_likelihood(method::YeoJohnson, 饾惐, 位)
    N = length(饾惐)
    饾惒 = transform(method, float.(饾惐), 位)
    蟽虏 = var(饾惒, corrected = false)
    c = sum(sign.(饾惐) .* log.(abs.(饾惐) .+ 1))
    llf = -N / 2.0 * log(蟽虏) + (位 - 1) * c
    @debug "YeoJohnson: 位 = $位 => 蟽虏=$蟽虏, c=$c, llf=$llf"
    return llf
end
